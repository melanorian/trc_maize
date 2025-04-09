# Goal:   This filters raw count tables extracted from dds based on choosen 
#         filtering criteria
# Input:  dds, choose filtering criteria
# Output: DDS_filtered, CTSNorm_filtered.csv, CPM_filtered.csv, vst_filtered

# Note: This script involves Exploratory elements and is intended to be run section
#       by section where choices or input is required you will encounter CAPITAL 
#       letters

# clean environment
rm(list = ls())

while (dev.cur() != 1) {
  dev.off()}

# SET FILTER CRITERIA
cutoffs <- c(1, 5, 10)
percentages <- c(0.3, 0.4, 0.5, 0.6, 0.7)

# 1. base set up install/load libraries, setwd, read files ----
# 1.1 load libraries
library(DESeq2)
library(reshape2)
library(tidyverse)
library(purrr)

# 1.2 SET working directory/ Prefixes/ Dynamic dir_out
# 1.2 Naming
# 1.1 load def_variables.R and def_functions.R
source("/media/wolf2/TKI2019/melanie/Code-transcriptomics-analysis/Pipeline/def_variables.R")
source("/media/wolf2/TKI2019/melanie/Code-transcriptomics-analysis/Pipeline/def_functions.R")

# 1.2 Set wd
setwd(dir_base)

# 1.3 Def in/out direcory

# 1.3.1 Load DDS
dir_in <- paste0(dir_base, "DDS/")

base_name <- f_make_name(select_genes = select_genes,
                         remove_samples = remove_samples)

name_dds <- paste("DDS", base_name, sep = "_")
path_dds <- paste0(dir_in, name_dds, ".RData")
dds <- f_assign_load(path_dds)

# 1.3.2 Dir out
dir_out_base <- paste0(dir_base, "Input_filter/")
dir_out <- paste0(dir_out_base, base_name)


if (!file.exists(dir_out_base)) {
  dir.create(dir_out_base)
}

if (!file.exists(dir_out)) {
  dir.create(dir_out)
}

# 2. Extract counts/metadata from dds
cts <- counts(dds)
coldata <- colData(dds) # METADATA

rownames(coldata) <- f_rename(rownames(coldata))
colnames(cts) <- f_rename(colnames(cts))

# 3. Filter the input cts using selected criteria
# 3.1 Function to filter cts based on overall rowsum

f_filter_counts <- function(cts, cutoff_count, cutoff_percentage) {
  # Select for only P1 samples by substring
  samples_oi <- grep("^t", colnames(cts))
  # Create a logical cts matrix for the samples of interest indicating if they are above the cutoff
  logical_cts <- cts[, samples_oi] > cutoff_count
  # Calculate how often soi exceed cutoff
  count_true <- rowSums(logical_cts)
  # Calculate threshold to fit filter criteria
  threshold <- round(length(samples_oi) * cutoff_percentage) # round number as we don't work with decimals
  # Extract list of relevant rownames
  filtered_cts <- cts[count_true >= threshold, ]
  return(filtered_cts)
}


# 3.2 use 3.1 in Nested lapply to iterate over cutoff values and percentages

filtered_list <- lapply(cutoffs, function(count_cutoff) {
  
  lapply(percentages, function(percentage_cutoff) {
    
    # Apply the function with the current combination of cutoff values
    filtered_cts <- f_filter_counts(cts = cts, 
                                    cutoff_count = count_cutoff, 
                                    cutoff_percentage = percentage_cutoff
                                    )
    
    # Return the filtered matrix and name it
    list(name = paste0("count_", 
                       count_cutoff, 
                       "_perc_", 
                       percentage_cutoff
                       ), 
         cts = filtered_cts
         )
  })
})

# 3.3 Flatten list and name the elements
flat_filter_list <- do.call(c, filtered_list)
names(flat_filter_list) <- sapply(flat_filter_list, function(x) x$name)

# 3.4 Remove the name element inside each list to clean up
flat_filter_list <- lapply(flat_filter_list, function(x) x$cts)

# reassigne to cts_filt for simplicity
cts_filt <- flat_filter_list
rm(filtered_list, flat_filter_list)

# 3.5 Print the filtered counts list
cat(names(cts_filt))

# 4. Generate new input files for downstream analysis with filtered output
# 4.1 calculate new CPMs
# Ensure each element in cts_filt is a matrix
cts_filt_mtx <- lapply(cts_filt, function(x) as.matrix(x))
CPM_filt_list <- lapply(cts_filt_mtx, f_CPM)

# 5.Re-construct dds
# 5.1 Prepare for constructing DESeq2 object ----
# 5.2 Define factor of interest
coldata <-  f_define_factor(coldata, factor_experiment[1])
#coldata <-  f_define_factor(coldata, factor_experiment[2])

# 5.3 Ensure factor_experiment is character
factor_experiment <- as.character(factor_experiment) 

# 5.6 If there is more than one factor...
# Construct the formula dynamically
if (length(factor_experiment) == 1) {
  # For a single factor, use 3.4 approach
  design_formula <- as.formula(paste("~", factor_experiment[1]))
} else {
  # If there is more than one factor, use 3.5 approach
  additional_terms <- factor_experiment[2]
  interaction_term <- paste(factor_experiment, collapse = ":")
  design_formula <- as.formula(paste("~", 
                                     factor_experiment[1],
                                     "+", 
                                     additional_terms,
                                     "+",
                                     interaction_term))
}

print(design_formula)

# 5.7 Function to generate dds dynamically
f_dds <- function(cts, coldata, design){
  
  dds_filterd <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = coldata,
    design = design_formula
  )
  return(dds_filterd)
}

# 5.8 Generate a list of dds 
dds_list <- lapply(cts_filt,
                   f_dds,
                   coldata = coldata,
                   design = design_formula
                   )

cat("generated filtered DDS...")

# 6. Extract normalised counts ----
# 6.1 Calculate size factors for extracting DESeq normalised counts
dds_list <- lapply(dds_list, estimateSizeFactors)

cat("calculated filtered size factors...")

# 6.2 Extract normalised counts
cts_den_norm <- lapply(dds_list, counts, normalized = TRUE) 

cat("generated DESeq2 normalised counts...")

# 6.3 Extract vst
cts_vst_list <- lapply(dds_list, vst) # default blind, true 
# 
# # Apply vst function to each DESeqDataSet object in the list
# cts_vst_list <- lapply(seq_along(dds_list), function(i) {
#   # Print message indicating the current DESeqDataSet object being processed
#   cat("Processing DESeqDataSet object at index:", i, "\n")
#   
#   # Apply vst function
#   vst(dds_list[[i]])
# })

vst_list <- lapply(cts_vst_list, assay)

### SAVE 
name_list <- names(cts_filt)

# Save filtered CPM counts
map2(.x = cts_filt,
     .y = paste0("CTS_", name_list),
     .f = f_safe_df,
     dir = dir_out,
     gois = "",
     part = "P1B4")

# Save filtered CPM counts
map2(.x = CPM_filt_list,
     .y = paste0("CPM_", name_list),
     .f = f_safe_df,
     dir = dir_out,
     gois = "",
     part = "P1B4")

# Function to save .RData = dds
f_save_rdata <- function(dds, name_dds, dir_save){
  save(dds,
       file = paste0(dir_save, 
                     "/",
                     "DDS_filtered_",
                     name_dds, 
                     ".RData")
       )
  
  paste("dds saved as:", paste0(dir_save, 
                               "/",
                               "DDS_filtered_",
                               name_dds, 
                               ".RData")
        )
}

# Save dds by calling functions
map2(.x = dds_list,
     .y = name_list,
     .f = f_save_rdata,
     dir_save = dir_out)


# Save DESEq2 normalised counts
map2(.x = cts_den_norm,
     .y = paste0("DESeqNorm", name_list),
     .f = f_safe_df,
     dir = dir_out,
     gois = "",
     part = "P1B4")
 
# Save DESEq2 normalised counts
map2(.x = vst_list,
     .y = paste0("VST_", name_list),
     .f = f_safe_df,
     dir = dir_out,
     gois = "",
     part = "P1B4")


# #################### Plott count distribution ###############################
# 
# # 2.2 Melt df into long format for plotting
# melted_cts <- melt(cts_filt$count_10_perc_0.7, 
#                    id.vars = colnames(cts_filt$count_10_perc_0.7), 
#                    variable.name = rownames(cts_filt$count_10_perc_0.7))
# 
# colnames(melted_cts) <- c("gene", "time", "expression")
# 
# # 2.3 Draw boxplot of log2(geneexpression/sample)
# melted_cts$log2 <- log2(melted_cts$expression +1) # due to large range of expression
# 
# p1 <- ggplot(melted_cts, aes(x = time, y = log2)) +
#   geom_boxplot() +
#   geom_point(position = position_jitter(width = 0.01), color = "grey") +
#   labs(x = "time", y = "log2(expression + 1)") +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# 
# # 2.4 Save the heatmap
# f_save_plot(dir = dir_out_plots, 
#             name = "/boxplot-rawreads",
#             filetype = "png",
#             plot = p1,
#             width = 10,
#             hight = 6)
# 
# ggsave(paste0(dir_out_plots, "boxplot-rawreads", ".png"), 
#        plot = p1, 
#        device = "png",
#        width = 8,
#        height = 5)
