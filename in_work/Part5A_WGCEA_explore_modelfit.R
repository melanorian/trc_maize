# Goal:   Identify clusters of co-expressed genes in input 
# Input:  
# Output: 

# 1.2 Load packages
library(DESeq2)
library(WGCNA)
library(tidyverse)
library(reshape)
library(ggplot2)
library(gridExtra)
library(parallel)

# PARAMETERS DEFINED IN def_variables.R

# 1.2 Define in and out directories
dir_in_dds <- paste0(dir_base, "DDS/")

dds_in_file <- paste0(
  dir_in_dds,
  "DDS_", 
  f_make_name(select_genes = select_genes, 
              remove_samples = remove_samples), "_",
  "DesNA",
  ".RData")

# out Directory
dir_WGCEA <- paste0(dir_base, "WGCEA/")

dir_out <- paste0(dir_WGCEA,
                  f_make_name(select_genes = select_genes,
                              remove_samples = remove_samples),
                  "/"
)

# 1.3 Make sure out directories exist
if (!file.exists(dir_WGCEA)) {
  dir.create(dir_WGCEA)
}

if (!file.exists(dir_out)) {
  dir.create(dir_out)
}


# 2. Load input data as dds
# 2.1 Check if dds file exists
if(file.exists(dds_in_file)){
  cat("Input DDS for WGCNA found as",
      dir_in_dds, "\n", 
      "Continue...")
  
} else {
  cat("Input DDS for WGCNA NOT! found as specified:",
      dir_in_dds, "\n", 
      "STOP...")
}

# 2.2 Load dds file 
load(dds_in_file)

# 2.3 extract coldata 
coldata <- colData(dds)

# 3. Filter input data based on min count in percentage of samples ----
# Note: already filtered for DGE but, apply stricter filter

# 3.1 Function to filter dds 
f_filter_dds <- function(dds, min_count, min_perc) {
  # Keep genes that have at least `min_count` in at least `min_samples` samples
  min_samples <- round(nrow(coldata)*min_perc)
  keep <- rowSums(counts(dds) >= min_count) >= min_samples
  dds_filtered <- dds[keep, ]
  
  return(dds_filtered)
}

# 3.2 Apply filter function
# for (count in count_range) { # activate for loop to run over range of counts
#   for (percent in percent_range) {

# file <- "_VST_dpi1_count"

dds_filtered <- f_filter_dds(
  dds = dds,
  min_count = filter_min_count,
  min_perc = filter_min_percent_samples
)

# 3.3 Extract vst data for WGCNA 
vst <- assay(vst(dds_filtered))

cat("Extracted filterd VST remaining genes:", nrow(vst), "/", nrow(dds))

# 3.4 Assign matrix to variable name "den"
den <- vst



# 4. Run WGCEA
# 4.1 Transform input data to match WGCEA requirements
#     WGCEA columns "gene probes", rows "treatments"
input_mtx <- t(den) 

# 4.2 Define nr threads, depends on available threads
allowWGCNAThreads(nThreads = round(detectCores()*0.2)) 

# 4.3 Define range of power to explore
powers <- c(c(1:10), seq(from = 12, to = 25, by = 2))

# 4.4 Call the network topology analysis function - takes a while
sft <- pickSoftThreshold(
  input_mtx,             # <= Input data
  blockSize = BLOCKSIZE,
  powerVector = powers,
  verbose = 5
)

# 4.5 Visualise the results & save to explore optimal soft power threshold
sft_df <- data.frame(sft$fitIndices) %>%
  dplyr::mutate(model_fit = -sign(slope) * SFT.R.sq)

f_plot_power <- function(sft_df) {
  ggplot(sft_df, aes(x = Power, y = model_fit, label = Power)) +
    geom_point() +
    geom_text(nudge_y = 0.1) +
    # WGCNA recommends R^2 cutoff
    geom_hline(yintercept = 0.80, col = "red") +
    # Just in case our values are low, we want to make sure we can still see the 0.80 level
    ylim(c(min(sft_df$model_fit), 1.05)) +
    xlab("Soft Threshold (power)") +
    ylab("Scale Free Topology Model Fit, signed R^2") +
    ggtitle("Scale independence") +
    # This adds some nicer aesthetics to our plot
    theme_classic()
}

f_plot_connection <- function(sft_df){
  ggplot(sft_df, aes(x = Power, y = mean.k.)) +
  geom_line(color = "grey") +  # Draw the mean connectivity line
  geom_point(color = "grey") +  # Add points on the line
  geom_text(aes(label = round(mean.k., 1)), 
            nudge_y = 0.1, 
            color = "red") +  # Label points with mean connectivity values
  labs(x = "Soft Threshold (power)", 
       y = "Mean Connectivity", 
       title = "Mean Connectivity") +
    theme_classic()
}

plot_power <- f_plot_power(sft_df = sft_df)
plot_connect <- f_plot_connection(sft_df = sft_df)
grid_power_connect <- arrangeGrob(plot_power, plot_connect, ncol = 2)

### AD TO NAME COUNT AND PERC
ggsave(
  filename = paste0(
    dir_out, 
    filter_min_count,
    filter_min_percent_samples,
    "grid_power",
    ".png"),
  plot =  grid_power_connect
  )

cat("Plotted optimal soft power thrashold for WGCNA at:", "\n",
    paste0(
      dir_out, 
      filter_min_count,
      filter_min_percent_samples,
      "grid_power",
      ".png"))
