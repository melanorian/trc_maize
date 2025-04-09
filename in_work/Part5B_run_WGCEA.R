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

cat("Extracted filterd VST remaining genes:", nrow(vst), "/", nrow(dds), "\n")

# 3.4 Assign matrix to variable name "den"
den <- vst



# 4. Run WGCEA
# 4.1 Transform input data to match WGCEA requirements
#     WGCEA columns "gene probes", rows "treatments"
input_mtx <- t(den) 

# 4.2 Define nr threads, depends on available threads
allowWGCNAThreads(nThreads = round(detectCores()*0.2)) 

# # 4.3 Define range of power to explore
# powers <- c(c(1:10), seq(from = 12, to = 25, by = 2))
# 
# # 4.4 Call the network topology analysis function - takes a while
# sft <- pickSoftThreshold(
#   input_mtx,             # <= Input data
#   blockSize = BLOCKSIZE,
#   powerVector = powers,
#   verbose = 5
# )
# 
# # 4.5 Visualise the results & save to explore optimal soft power threshold
# sft_df <- data.frame(sft_test$fitIndices) %>%
#   dplyr::mutate(model_fit = -sign(slope) * SFT.R.sq)
# 
# f_plot_power <- function(sft_df) {
#   ggplot(sft_df, aes(x = Power, y = model_fit, label = Power)) +
#     geom_point() +
#     geom_text(nudge_y = 0.1) +
#     # WGCNA recommends R^2 cutoff
#     geom_hline(yintercept = 0.80, col = "red") +
#     # Just in case our values are low, we want to make sure we can still see the 0.80 level
#     ylim(c(min(sft_df$model_fit), 1.05)) +
#     xlab("Soft Threshold (power)") +
#     ylab("Scale Free Topology Model Fit, signed R^2") +
#     ggtitle("Scale independence") +
#     # This adds some nicer aesthetics to our plot
#     theme_classic()
# }
# 
# f_plot_connection <- function(sft_df){
#   ggplot(sft_df, aes(x = Power, y = mean.k.)) +
#   geom_line(color = "grey") +  # Draw the mean connectivity line
#   geom_point(color = "grey") +  # Add points on the line
#   geom_text(aes(label = round(mean.k., 1)), 
#             nudge_y = 0.1, 
#             color = "red") +  # Label points with mean connectivity values
#   labs(x = "Soft Threshold (power)", 
#        y = "Mean Connectivity", 
#        title = "Mean Connectivity") +
#     theme_classic()
# }
# 
# plot_power <- f_choose_power(sft_df = sft_df)
# plot_connect <- f_plot_connection(sft_df = sft_df)
# grid_power_connect <- arrangeGrob(plot_power, plot_connect, ncol = 2)
# 
# ggsave(
#   filename = paste0(
#     dir_out, 
#     "grid_power",
#     ".png"),
#   plot =  grid_power_connect
  # )

# 4. Run WGCNA using the power previously defined
f_run_block <- function(input_matrix, power, dir_log) {
  # Force it to use WGCNA cor function (fix a namespace conflict issue)
  temp_cor <- cor
  cor <- WGCNA::cor
  
  # Open a connection to a log file
  sink(paste0(dir_log, "logfile.txt"), append = TRUE)
  
  # log which parameters
  cat("input matrix:", deparse(substitute(input_matrix)), "\n")
  cat("power:", power, "\n")
  cat("Network construction blockwiseModules started..\n")
  
  # Run blockwiseModules function
  ntwk <- blockwiseModules(input_matrix,
                            power = power, # soft threshold for network construction
                          networkType = "signed", # Important parameter!
                            deepSplit = 2,
                            pamRespectsDendro = FALSE,
                            minModuleSize = 30,
                            maxBlockSize = BLOCKSIZE, # What size chunks (how many genes) the calculations should be run in
                            reassignThreshold = 0,
                            mergeCutHeight = 0.25,
                            #randomSeed = 1234, # there's some randomness associated with this calculation
                            # so we should set a seed
                            TOMType = "signed", # topological overlap matrix
                            saveTOMs = TRUE,
                            saveTOMFileBase = "ER",
                            numericLabels = TRUE,  #Let's use numbers instead of colors for module labels
                            verbose = 3,
                            )
  
  cat("Network construction finished...")
  
  # Close the connection to the log file
  sink()
  
  # Restore the original cor function
  cor <- temp_cor
  
  return(ntwk)
}

ntwk <- f_run_block(input_matrix = input_mtx,
                    power = POWER,
                    dir_log = dir_out)

readr::write_rds(
  ntwk,
  file = paste0(
    dir_out, 
    paste0(
      "WGCNA_RES_", 
      "count", filter_min_count, 
      "perc", filter_min_percent_samples, 
      "_power", POWER, 
      "_block", BLOCKSIZE, 
      ".RDS"
    )
  )
)


# 5. Function to plot results of WGCEA as dendogram
f_plot_dendo <- function(ntwk, dir_save, filetype) {
  # Convert labels to colors for plotting
  mergedColors <- labels2colors(ntwk$colors)
  
  # Open a PNG graphics device
  file_path <- paste0(
    dir_save, 
    paste0(
      "WGCNA_RES_", 
      "count", filter_min_count, 
      "perc", filter_min_percent_samples, 
      "_power", POWER, 
      "_block", BLOCKSIZE,
      "DENDO_Cluster", 
      ".", 
      filetype))
  
  png(file_path, width = 10, height = 5, units = "in", res = 300)
  
  # Plot the dendrogram and the module colors underneath
  plotDendroAndColors(
    ntwk$dendrograms[[1]],
    mergedColors[ntwk$blockGenes[[1]]],
    "Module colors",
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05
  )
  
  # Close the graphics device
  dev.off()
}

f_plot_dendo(ntwk = ntwk,
             dir_save = dir_out,
             filetype = "png")

# 6. Function to make a df where each gene is sorted into a Module
f_modules_df <- function(ntwk) {
  # Convert labels to colors for creating the data frame
  moduleColors <- labels2colors(ntwk$colors)
  
  # Create a data frame with gene IDs and their corresponding module colors
  module_df <- data.frame(
    gene_id = names(ntwk$colors),
    colors = moduleColors
  )
  
  # Return the module data frame
  return(module_df)
}


module_df <- f_modules_df(ntwk = ntwk)

# 6.2 Safe in a table
write.table(module_df,
           file = paste0(
             dir_out, 
              paste0(
                "WGCNA_RES_", 
                "count", filter_min_count, 
                "perc", filter_min_percent_samples, 
                "_power", POWER, 
                "_block", BLOCKSIZE,
                "GENE_MODULES", 
                ".", 
                "txt")),
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

cat("module dataframe safed at:", dir_out)


############################ CONTINUE HERE ##############################
############################ CONTINUE HERE ##############################

# 7. Function to extract Eigengenes for each cluster and sample order based on
# 7. similarity
#ntwk <- ntwk
f_module_eigengenes <- function(input_mtx, ntwk) {
  mergedColors <- labels2colors(ntwk$colors)
  
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

# 7.2 Apply function
eigen <- f_module_eigengenes(input_mtx = input_mtx, ntwk = ntwk)

# 8. Make a plot to analyse Module-trait Relationships
f_module_rel_plot <- function(MEs0) {
  # Tidy data
  mME <- MEs0 %>%
    pivot_longer(-treatment) %>%
    mutate(
      name = gsub("ME", "", name),
      name = factor(name, levels = colnames(MEs0))
    )
  
  # Plot data
  plot <- ggplot(mME, aes(x = treatment, y = name, fill = value)) +
    geom_tile() +
    theme_bw() +
    scale_fill_gradient2(
      low = "blue",
      high = "red",
      mid = "white",
      midpoint = 0,
      limit = c(-1, 1)
    ) +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(title = "Module-trait Relationships", y = "Modules", fill = "corr")
  
  return(plot)
}

p_mod_traitrel <- f_module_rel_plot(MEs0 = eigen)

f_save_plot(dir = dir_out,
            name = paste0(
              "ModuleTraitRel",
              "count", filter_min_count, 
              "perc", filter_min_percent_samples, 
              "_power", POWER, 
              "_block", BLOCKSIZE), 
            filetype = "png",
            p_mod_traitrel,
            width = 8,
            hight = 5
            )
#   }  ######### switch on for for loop
# }

# ________________________________________________________________________-
# 
# # # 9. Plot the expression patterns of selected gene-subsets
f_module_expression <- function(modules, module_df, expression_mtx){

  submod <- module_df %>%
    subset(colors %in% modules)

  row.names(module_df) = module_df$gene_id

  # subset expression matrix for modules of interest
  subexpr <- expression_mtx[submod$gene_id,]

  return(subexpr)
}

f_plot_module <- function(subexpr, module_df, y_axis_label) {

  # Set row names before subsetting
  row.names(module_df) <- module_df$gene_id

  # Convert the expression data to a long format for ggplot
  submod_df <- data.frame(subexpr) %>%
    mutate(
      gene_id = row.names(.)
    ) %>%
    pivot_longer(-gene_id) %>%
    mutate(
      module = module_df[gene_id,]$colors
    )

  # Plot expression data
  p <- submod_df %>% ggplot(aes(x = name, y = value, group = gene_id)) +
    geom_line(aes(color = module), alpha = 0.2) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) +
    facet_grid(rows = vars(module)) +
    labs(x = "treatment", y = y_axis_label)

  return(p)
}

col_oi <- c("grey", "yellow", "black", "brown", "pink", "magenta", "red", "turquoise", "green", "blue", "purple")

for (colour in col_oi){

subexpr <- f_module_expression(modules = colour,
                               module_df = module_df,
                               expression_mtx = den)

p_expr_module <- f_plot_module(subexpr = subexpr,
                               module_df = module_df,
                               y_axis_label = "vst expression")


f_save_plot(dir = dir_out,
            name = paste0("expr_module", paste0(
              "count", filter_min_count, 
              "perc", filter_min_percent_samples, 
              "_power", POWER, 
              "_block", BLOCKSIZE), colour),
            filetype = "png",
            plot = p_expr_module,
            width = 8,
            hight = 4)

 cat("module expression plot saved at:", dir_out)

}


# ___________________________________________________________________
# # Calculate z-scores for each gene's expression values
# subexpr_z <- t(scale(t(subexpr)))  # Transpose, scale, and transpose back
# subexpr_log <- t(log2(t(subexpr)))
# netwk

#______________________________________________

# # Export for e.g. cytoscape
# genes_of_interest = module_df %>%
#   subset(colors %in% modules_of_interest)
# 
# expr_of_interest = den[genes_of_interest$gene_id,]
# 
# # Only recalculate TOM for modules of interest 
# TOM = TOMsimilarityFromExpr(t(expr_of_interest),
#                             power = POWER)
# 
# # Add gene names to row and columns
# row.names(TOM) = row.names(expr_of_interest)
# colnames(TOM) = row.names(expr_of_interest)
# 
# edge_list = data.frame(TOM) %>%
#   mutate(
#     gene1 = row.names(.)
#   ) %>%
#   pivot_longer(-gene1) %>%
#   dplyr::rename(gene2 = name, correlation = value) %>%
#   unique() %>%
#   subset(!(gene1==gene2)) %>%
#   mutate(
#     module1 = module_df[gene1,]$colors,
#     module2 = module_df[gene2,]$colors
#   )
# 
# head(edge_list)
# 
# write.table(edge_list,
#             file = paste0(dir_out, "edgelist.tsv", filter),
#             row.names = TRUE,
#             col.names = TRUE,
#             sep = "\t")
# 
# cat("edge list saved at:", dir_out)
