# Goal:   This script explores the effect of filtering genes with low read counts
#         on the number of genes remaining in the dataset
# Input:  Requires a dds with samples and genes of interest
# Output: Heatmaps overlayed with number of genes remaining in the dataset 
#         appling defined combinations of filtering criteria

# 1. base set up install/load libraries, setwd, read files ----
# 1.1 load libraries
library(DESeq2)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(pheatmap)
library(tidyr)

dir_in <- paste0(
  dir_base, "Input_Matrix/",
  f_make_name(
    select_genes = select_genes,
    remove_samples = remove_samples),
  "/"
)

# Load CTS
f_load_if_absent_csv(
  dir_in = dir_in,
  name_env_variable = "coldata",
  match_str = "METADATA")

f_load_if_absent_csv(
  dir_in = dir_in,
  name_env_variable = "cts",
  match_str = "CTS")

###########################
# 1.3 Def in/out direcory

# 1.3.1 Load DDS
# dir_in <- paste0(dir_base, "DDS/")

base_name <- f_make_name(select_genes = select_genes,
                             remove_samples = remove_samples)

# name_dds <- paste("DDS", base_name, sep = "_")
# path_dds <- paste0(dir_in, name_dds, ".RData")
# dds <- f_assign_load(path_dds)

# 1.3.2 Dir out
dir_out_plots <- paste0(
  dir_base, 
  "EDA_basic", "/", 
  base_name
  )  

######
# ADD CHECK IF DIR EXIST
######

# # 2. Extract counts/metadata from dds
# cts <- counts(dds)
# coldata <- colData(dds) # METADATA

# 2.1 Define range of cutoff options to
cutoffs <- c(5:17)
percentages <- seq(0.1, 0.9, by = 0.1) # fraction of samples to meet count cutoff

exp_factor <- c(1:8) # timepoints

# 2.2 Melt df into long format for plotting
# melted_cts <- melt(cts, id.vars = colnames(cts), variable.name = rownames(cts))
# 
# colnames(melted_cts) <- c("gene", "time", "expression")

melted_cts <- as.data.frame(cts) %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    cols = -gene,            # Keep 'gene' column as id
    names_to = "sample",     # Extract sample names
    values_to = "expression" # Values become 'expression' column
  )

# Large range of expression, log transform
melted_cts$log2 <- log2(melted_cts$expression +1) # due to large range of expression

# 2.3 Draw boxplot of log2(geneexpression/sample)
p1 <- ggplot(melted_cts, 
             aes(x = sample, y = log2)
             ) +
  
  geom_boxplot(outlier.colour = "red", 
               outlier.shape = 16,
               outlier.size = 1
               ) +
  
  geom_point(position = position_jitter(
    width = 0.01),
    color = "grey",
    size = 1,
    alpha = 0.5
    ) + # Transparent jitter points
  
  labs(x = "Sample", 
       y = "log2(expression + 1)"
       ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45, 
      hjust = 1)
    )
  

# 2.4 Save the heatmap
f_save_plot(dir = dir_out_plots, 
            name = "/boxplot-rawreads",
            filetype = "png",
            plot = p1,
            width = 15,
            hight = 6)


# 3. Plot number of genes above a certain threshold
# 3.1 Convert count matrix into df
df_cts <- as.data.frame(cts)

# 3.2 list samples 
samples <- colnames(df_cts)

# 3.3 Define function to count number of times genes are expressed above a cutoff
f_count_nrgenes <- function(cts_df, sample, cutoff){
  column <- cts_df[[sample]]
  instances <- sum(column > cutoff)
  return(instances)
}

# 3.5 Create empty df to store number of times genes are expressed above a cutoff
counts_nrgenes <- data.frame(matrix(NA, nrow = length(samples), ncol = length(cutoffs)))
rownames(counts_nrgenes) <- samples
colnames(counts_nrgenes) <- cutoffs

# 3.6 Run function f_count_nrgenes over all defined count cutoffs 
for (cutoff in cutoffs) {
  # Apply the function to each sample
  counts <- sapply(samples, function(sample) f_count_nrgenes(df_cts, sample, cutoff))
  # Store the results in the dataframe
  counts_nrgenes[, as.character(cutoff)] <- counts
}

# 3.7 Make uniform colnames for plotting data & transform to long format
counts_nrgenes$Sample <- rownames(counts_nrgenes)

melted_counts <- melt(counts_nrgenes, 
                      id.vars = "Sample", 
                      variable.name = "Cutoff", 
                      value.name = "Count"
                      )

f_safe_df(df = counts_nrgenes,
          name_df = "nr_genes_above_cutoff",
          dir = dir_out_plots,
          gois = "all",
          part = "P1B3")

# 3.8 Plot by timepoint --> ensure that all samples of the same timepoint have the same name
substring_names <- sub("_(?!.*_).*", "", melted_counts$Sample, perl = TRUE)
melted_counts$Sample <- substring_names

# 3.9 for using split function cutoffs cannot be numeric, so add the word cutoff
melted_counts$Cutoff <- as.character(paste0(
  melted_counts$Cutoff,
  "_cutoff" 
  )
  )

head(melted_counts)
# 3.10 Split df into individual df for plotting, split by cutoff
melted_counts_list <- split(melted_counts, melted_counts$Cutoff)

# 3.11 Plotting, the real deal
# 3.11.1 Define function to generate plot & safe output

f_nrcounts <- function(df, name, variable) {
  
  # Draw heatmap
  p2  <- # Assuming your dataframe is named df
    ggplot(df, 
           aes(x = Sample, y = !!as.name(variable))
           ) +
    geom_boxplot() +
    geom_jitter(width = 0.1, 
                height = 0, 
                color = "grey", 
                alpha = 0.5
                ) +
    geom_hline(yintercept = 25000, 
               linetype = "dashed", 
               color = "black"
             ) +
    geom_hline(yintercept = 22500, 
               linetype = "dashed", 
               color = "black"
    ) +
    geom_hline(yintercept = 20000, 
               linetype = "dashed", 
               color = "black"
    ) +
    ggtitle(name) +
    labs(x = "Sample",
         y = paste("nr genes >", name)
         ) +
    coord_cartesian(ylim = c(15000, 30000)) +
    theme_classic() +
    theme(
      text = element_text(size = 15),
      axis.text.x = element_text(angle = 45, hjust = 1))

  return(p2)
}



# 3.11.2 Extract names of matrices in subsetted_mtx_list
melted_names <- names(melted_counts_list)

# 3.11.3 Apply the function to each matrix in subsetted_mtx_list
p2_list <- map2(
  .x = melted_counts_list,
.y = melted_names,
  .f = ~f_nrcounts(
    df = .x,
    name = .y,
    variable = "Count"
  )
)


p2_list<- p2_list[order(as.numeric(
  gsub("_cutoff", "", names(p2_list))))]

# Optionally, save the arranged grid plot as a single image
p2_combined <- arrangeGrob(
  grobs = p2_list,
  ncol = round(sqrt(length(p2_list)))
  )

ggsave(filename = paste0(dir_out_plots, 
                         "/boxplot_grid_cutoff",
                         paste0(
                           "_", 
                           as.character(min(cutoffs)), 
                           "to",
                           as.character(max(cutoffs))
                             ),
                         ".png"), 
       plot = p2_combined, 
       width = 5*round(sqrt(length(p2_list))), 
       height =6*round(sqrt(length(p2_list)))
       )


# Save the heatmap
ggsave(paste0(dir_out_plots, 
              "/", 
              "plot", 
              name, 
              ".png"
), 
plot = p2, 
device = "png",
width = 5,
height = 3.5
)

# 4. Filter  ----
# 4.1 Filter by simple rowsum (raw counts)
# 4.1.1 Function for filtering by rowsum
f_filter_rowsum <- function(count_table, rowsum_cutoff){
  row_sum <- rowSums(count_table)
  filtered_cts <- count_table[row_sum >= rowsum_cutoff, ]
  return(filtered_cts)
}

# 4.1.2 Run function for filtering by rowsum and name adequately
cts_filtered_list <- lapply(filter, 
                            f_filter_rowsum, 
                            count_table = cts
                            )

cts_filtered_names <- paste(paste("CTSnorm", 
                                  subset_out, 
                                  "filtered", 
                                  sep = "_"
                                  ), 
                            paste("rowsum", 
                                  filter, 
                                  sep = "_"
                                  )
                            )

names(cts_filtered_list) <- cts_filtered_names

# 4.2 Explore filtering for genes with >X counts in at least Y% of samples at Zexp_factor ----

# 4.2.1 Define range of interest for cutoffs & range timepoints
new_cols <- f_rename(colnames(cts))
colnames(cts) <- new_cols

# 4.2.2 function to extract column index with samples corresponding to a specific timepoint
f_colindex_exp_factor <- function(cts, exp_factor) {
  colnames <- colnames(cts)
  colindex_exp_factor <- grep(paste0("^t", exp_factor), colnames)
  print(exp_factor)
  print(colindex_exp_factor)
  return(colindex_exp_factor)
}

# 4.2.3 Define a function which allows to filter the cts for defined cutoffs 
f_filter_counts <- function(cts, samples_oi, cutoff_count, cutoff_percentage) {
  # Create a logical cts matrix for the samples of interest indicating if they are above the cutoff
  logical_cts <- cts[, samples_oi] > cutoff_count
  # Calculate how often soi exceed cutoff
  count_true <- rowSums(logical_cts)
  # Calculate threshold to fit filter criteria
  threshold <- round(length(samples_oi) * cutoff_percentage) # round number as we don't work with decimals
  # Extract list of relevant rownames
  filtered_genes <- rownames(cts[count_true >= threshold, ])
 
  return(filtered_genes)
}

# 4.2.4 Initialize an empty list to store explor_cutoffs_list
explor_cutoffs_list <- list()

# 4.2.5 Nested loop to run all comparison for each exp_factor value
for (current_exp_factor in exp_factor) {
  # Get column indices for the current exp_factor
  samples_oi <- f_colindex_exp_factor(cts= cts, exp_factor = current_exp_factor)
  # Generate all possible combinations of cutoffs and percentages for the exp_factor (DEFINE range above)
  combinations <- expand.grid(cutoff_count = cutoffs, cutoff_percentage = percentages)
  # Initialize a list to store the explor_cutoffs_list for the current exp_factor
  exp_factor_results <- list()
  
  # Loop over each combination of cutoffs and percentages by going over rows in combinations
  for (i in 1:nrow(combinations)) {
    # Get the current combination
    current_combination <- combinations[i, ]
    # Apply the f_filter_counts function with the current combination
    filtered_counts <- f_filter_counts(cts = cts, 
                                       samples_oi = samples_oi,
                                       cutoff_count = current_combination$cutoff_count, 
                                       cutoff_percentage = current_combination$cutoff_percentage)
    # Store the filtered counts in the explor_cutoffs_list list
    exp_factor_results[[paste0("cutoff_", current_combination$cutoff_count, "_percentage_", current_combination$cutoff_percentage)]] <- filtered_counts
  }
  
  # Store the explor_cutoffs_list for the current exp_factor
  explor_cutoffs_list[[paste0("exp_factor_", current_exp_factor)]] <- exp_factor_results
}

# 4.2.6 Define function to process each exp_factor result into a matrix
f_process_exp_factor_results <- function(exp_factor_results) {
  # Get the lengths of each filtered result for a specific exp_factor
  lengths <- sapply(exp_factor_results, length)
  # Add a new column to combinations containing the lengths of results
  combinations$length_of_results <- lengths
  # Transform the results data frame to the desired format
  filter_stats <- combinations %>%
    pivot_wider(names_from = cutoff_count, values_from = length_of_results) %>%
    column_to_rownames(var = "cutoff_percentage")
  
  return(filter_stats)
}

# 4.2.7 Apply the function to each element of explor_cutoffs_list
explore_cutoffmtx_list <- lapply(explor_cutoffs_list, f_process_exp_factor_results)

# 5. draw heatmap of exploration for comprehensive insight https://slowkow.com/notes/pheatmap-tutorial/#uniform-breaks
# 5.1 Define colour scheme 
paletteLength <- 50
myColor <- colorRampPalette(c("#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF", "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026"))(paletteLength)

# 5.2 Define a function to generate and save heatmaps
f_heatmap <- function(mtx, names) {
  # # File name
  # filename <- paste(pre, pt, subset_out, "overview-heatmap", "exp_factor", names , sep = "_")
  
  # Draw heatmap
  h1 <- pheatmap(mtx,
                 display_numbers = T,
                 number_format = "%.0f",
                 number_color = "black",
                 color = myColor,
                 cluster_cols = F,
                 cluster_rows = F,
                 treeheight_row = 0,
                 angle_col = 0)
  
  # Save the heatmap
  ggsave(paste0(dir_out_plots, 
                "/", 
                names, 
                "_",
                "heatmap-filtering", 
                ".png"
                ), 
         plot = h1, 
         device = "png",
         width = 5,
         height = 3)
  
  # Close graphics device
  dev.off()
  
  return(h1$gtable)
}

# 5.3 Extract the name of the matrix elements saved in explore_cutoffmtx_list
matrix_names <- names(explore_cutoffmtx_list)

# 5.4 Apply the heatmap function to each matrix in explore_cutoffmtx_list, using the matrix names to acces each element
h1_list <- lapply(matrix_names, function(name) {
  f_heatmap(explore_cutoffmtx_list[[name]], name)
})

h1_combined <- arrangeGrob(grobs = h1_list, ncol = 2)
ggsave(filename = paste0(dir_out_plots, "/heatmap_cutoff.png"), 
       plot = h1_combined, 
       width = 10, 
       height = 12
)

# 6. at least >X reads in Y% samples at 1 exp_factor & at least >Z reads in A% samples all exp_factor -1 exp_factor

#############------------

# 7. Filter based on 1 exp_factor  

f_filterby_exp_factor_lowest <- function(cts, exp_factor_lowest, cutoff_count, cutoff_percentage) {
  exp_factor_name <- paste0("^t", exp_factor_lowest)
  samples_oi <- grep(exp_factor_name, colnames(cts))
  # Create a logical cts matrix for the samples of interest indicating if they are above the cutoff
  logical_cts <- cts[, samples_oi] > cutoff_count
  # Calculate how often soi exceed cutoff
  count_true <- rowSums(logical_cts)
  # Calculate threshold to fit filter criteria
  threshold <- round(length(samples_oi) * cutoff_percentage) # round number as we don't work with decimals
  # Extract list of relevant rownames
  filtered_genes <- rownames(cts[count_true >= threshold, ])
  
  return(filtered_genes)
}

trial <- f_filterby_exp_factor_lowest(cts = cts, 
                      exp_factor_lowest = "1", 
                      cutoff_count = 1,
                      cutoff_percentage = 0.5)

f_filter_exp_factor_others <- function(cts, exp_factor_lowest, cutoff_count, cutoff_percentage) {
  exp_factor_name <- paste0("^t", exp_factor_lowest)
  samples_oi <- grep(exp_factor_name, colnames(cts))
  # Create a logical cts matrix for the samples of interest indicating if they are above the cutoff
  logical_cts <- cts[, samples_oi] > cutoff_count
  # Calculate how often soi exceed cutoff
  count_true <- rowSums(logical_cts)
  # Calculate threshold to fit filter criteria
  threshold <- round(length(samples_oi) * cutoff_percentage) # round number as we don't work with decimals
  # Extract list of relevant rownames
  filtered_genes <- rownames(cts[count_true >= threshold, ])
  
  return(filtered_genes)
}


####----------------------

# 6. Filter cts based on exploration 
# 4.2.3.1 Adapt f_filter_counts function to return cts df 
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

# 4.2.3.2 Run function over cutoff/percentage combinations

mtx_rowfilter <- matrix(0, 
                        nrow = length(percentages), 
                        ncol = length(cutoffs),
                        dimnames = list(paste0(percentages), 
                                        paste0(cutoffs)))

# Using nested lapply to iterate over cutoff values and percentages
lapply(1:length(cutoffs), function(i) {
  lapply(1:length(percentages), function(j) {
    count_cutoff <- cutoffs[i]
    percentage_cutoff <- percentages[j]
    # Apply the function with the current combination of cutoff values
    filtered_cts <- f_filter_counts(cts = cts, cutoff_count = count_cutoff, cutoff_percentage = percentage_cutoff)
    # Store the number of rows in the result matrix
    mtx_rowfilter[j, i] <<- nrow(filtered_cts)
  })
})

print(mtx_rowfilter)


# Plot and safe results for easy interpretation
f_heatmap(mtx = mtx_rowfilter, names = "rowsum")
