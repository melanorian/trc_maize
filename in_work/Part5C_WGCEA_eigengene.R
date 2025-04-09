# Goal:   
# Input:  
# Output: 

# 1.2 Load packages
library(WGCNA)
library(DESeq2)
library(tidyverse)
library(reshape)
library(ggplot2)
library(gridExtra)

# PARAMETERS DEFINED IN def_variables.R

# 1.2 Define in and out directories
dir_in_dds <- paste0(dir_base, "DDS/")
dir_WGCEA <- paste0(dir_base, "WGCEA/")

dir_in_WGCNA <- paste0(dir_WGCEA,
                  f_make_name(select_genes = select_genes,
                              remove_samples = remove_samples),
                  "/"
)

dds_in_file <- paste0(
  dir_in_dds,
  "DDS_", 
  f_make_name(select_genes = select_genes, 
              remove_samples = remove_samples), "_",
  "DesNA",
  ".RData")

WGCNA_in_file <- paste0(dir_in_WGCNA,
                        paste0(
                          "WGCNA_RES_", 
                          "count", filter_min_count, 
                          "perc", filter_min_percent_samples, 
                          "_power", POWER, 
                          "_block", BLOCKSIZE, 
                          ".RDS"
                          )
)

# 1.3 Make sure out directories exist
if (!dir.exists(dir_in_dds)) {
  dir.create(dir_in_dds)
}

if (!dir.exists(dir_WGCEA)) {
  dir.create(dir_WGCEA)
}

if (!dir.exists(dir_in_WGCNA)) {
  dir.create(dir_in_WGCNA)
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

# 2.1 Check if WGCNA file exists
if(file.exists(WGCNA_in_file)){
  cat("Input WGCNA Results found as", "\n",
      WGCNA_in_file, "\n", 
      "Continue...")
  
} else {
  cat("WGCNA Results NOT! found as specified:", "\n",
      WGCNA_in_file, "\n", 
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

# 3.5 For WGCNA need to tranfrom col/row in matrix
input_mtx <- t(den)

# 4. Load WCGNA 
if (!exists("ntwk", envir = .GlobalEnv)) {
  ntwk <- readr::read_rds(WGCNA_in_file)
  cat("WCGNA loaded succesfully. Continue extracting Eigengenes...")
}else{
  cat("WCGNA already loaded. Continue extracting Eigengenes...")
}

# 5. Extract eigen modules from
# 5.1 Use functions from WGCENA 
f_module_eigengenes <- function(input_mtx, netwk) {
  mergedColors <- labels2colors(netwk$colors)
  # Get Module Eigengenes per cluster
  MEs0 <- moduleEigengenes(input_mtx, mergedColors)$eigengenes
  # Reorder modules so similar modules are next to each other
  MEs0 <- orderMEs(MEs0)
  # Generate module order without "ME" prefix
  module_order <- names(MEs0) %>% gsub("ME", "", .)
  colnames(MEs0) <- module_order
  # Add treatment names
  MEs0$treatment <- row.names(MEs0)
  return(MEs0)
}

# 5.2 eigengene funcion Apply function
eigen <- data.frame(f_module_eigengenes(input_mtx = input_mtx, netwk = ntwk))

# 6. Draw line plots for eigengenes
# 6.1 transform eigen to a long format 
eigen_long <- eigen %>%
  rownames_to_column(var = "sample") %>%
  gather(key = "color", value = "value", -sample)

# 6.2 Filter out rows with "treatment" in the color column
eigen_long <- eigen_long %>%
  filter(color != "treatment")

# 6.3 Separate the sample_group info from sample name
eigen_long <- eigen_long %>%
  mutate(sample_group = str_extract(sample, "T[0-9]+"))

# 6.4. Ensure the eigen values are numeric
eigen_long$value <- as.numeric(eigen_long$value)

# 6.5 Calculate median and standard deviation for each sample_group and color
summary_stats <- eigen_long %>%
  #filter(grepl("^t", sample)) %>%
  group_by(sample_group, color) %>%
  summarize(median_value = median(value), sd_value = sd(value))

# 6.6 Function to draw plots of eigengene values over time 
f_eigen_plot <- function(colour){
  # Filter statistics for the color "purple"
  colour_stats <- summary_stats %>%
  filter(color == colour)

  # Merge the colour statistics with the original data
  colour_data <- merge(eigen_long, colour_stats, by = "sample_group")
  
  e1 <- ggplot(colour_stats,
               aes(x = sample_group,
                   y = median_value)
               ) +
  geom_line(color = colour,
            size = 1.5, 
            group = 1
            ) +
  geom_point(color = colour, size = 3
             ) +
  geom_errorbar(aes(ymin = median_value - sd_value, 
                    ymax = median_value + sd_value), 
                width = 0.1
                ) +
  labs(title = paste("Eigengene module", colour), 
       x = "dpi", 
       y = "Median eigengene value"
       ) +
  coord_cartesian(ylim = c(-0.3, 0.3)
                  )+
  theme_minimal()

return(e1)
}

# 6.7 Run functio for eigen plots
modules <- colnames(eigen)
p_eigen <- lapply(modules, f_eigen_plot)

# 6.8 Combine all plots into a grid
grid_plot <- grid.arrange(grobs = p_eigen, ncol = round(sqrt(length(p_eigen))))

# 5. Safe Plots
# 5.1 Individual plots
for (i in seq_along(p_eigen)) {
  ggsave(paste0(dir_in_WGCNA,
                "EIGENGENES_",
                "count", filter_min_count, 
                "perc", filter_min_percent_samples, 
                modules[i], 
                ".png"
                ),
         plot = p_eigen[[i]]
         )
}

# 5.2 Save Grid
ggsave(paste0(dir_in_WGCNA, 
              "EIGENGENES_", 
              "count", filter_min_count, 
              "perc", filter_min_percent_samples, 
              "_power", POWER, 
              "_block", BLOCKSIZE,  
              "_GRID.png"
              ), 
       plot = grid_plot,
       width = 20,
       height = 15
     )

cat("Eigengene plot safed as:", "\n",
    paste0(dir_in_WGCNA, 
           "EIGENGENES_", 
           "count", filter_min_count, 
           "perc", filter_min_percent_samples, 
           "_power", POWER, 
           "_block", BLOCKSIZE,  
           "_GRID.png"
    )
)
