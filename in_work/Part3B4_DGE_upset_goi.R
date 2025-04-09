# Goal:   This script allows creating a binary matrix for DGE data in order to 
#         plot the results as an upsetplot, and identify which genes are uniquely expressed
#         at a specific timepoint or over a set number of timepoints
# Input:  results of DEG analysis in .RData format where res$contrast$direction 
# Output: Binary matrix, upsetplot, list of genes uniquely DE at a specific timepoint or 
#         expressed over all timepoints (possibly interesting candidates)

# 1. Set-up & import data ----
# Clean environment
rm(list = ls())

while (dev.cur() != 1) {
  dev.off()}

# 1.1 load libraries
library(dplyr)
library(purrr)
library(svglite)

# 1.1 load def_variables.R and def_functions.R
source("/media/wolf2/TKI2019/melanie/Code-transcriptomics-analysis/Pipeline/def_variables.R")
source("/media/wolf2/TKI2019/melanie/Code-transcriptomics-analysis/Pipeline/def_functions.R")

# 1.2 Set wd
setwd(dir_base)

# 1.3 define set name
set_name <- f_make_name(select_genes = select_genes,
                        remove_samples = remove_samples)


# 1.5 out directory
dir_goi <- paste0(dir_base, "DGE-goi", "/", 
                  set_name, "_", gois
)

name_sub <- paste(ref_gr, "filterLFC", 
                  paste0("FC", log2FC), 
                  paste0("PADJ", pval), 
                  sep = "_"
)

dir_in <- paste0(dir_goi, "/", name_sub)


# 1.5 Define out directory 
dir_out <- dir_in

# 2. Load files
# 2. Load binary matrix

binema_file <- f_list_files(dir = dir_in,
                            match_str = "_binema_",
                            file_type = "csv")

file_names <- gsub("P3B_|\\_.csv", "", binema_file)

mtx_binema <- lapply(binema_file, 
       f_read_csv,
       dir = paste0(dir_in, "/")
       )

cat("numbr of binemas_in:", length(mtx_binema))

names(mtx_binema) <- file_names
cat("names of binemas:", file_names)


# A. For biplot input: remove row names of each matrix in the list
biplot_list_in <- lapply(mtx_binema, function(x) {
  rownames(x) <- NULL
  return(x)
})

names(biplot_list_in) <- paste0(de_direction, "_biplot")


f_upset_plot <- function(biplot_in, direction) {
  name_safe_svg <- paste0(dir_out, "/", direction, "_upset.svg")
  
  # Create the UpSet plot
  upset_plot <- upset(data = biplot_in,
                      sets = rev(colnames(biplot_in)), # Customize based on desired order
                      nsets = length(colnames(biplot_in)),
                      keep.order = TRUE,
                      sets.x.label = "sum differentially expressed genes",
                      matrix.color = "#E69F00",
                      main.bar.color = "#E69F00",
                      mainbar.y.label = "number differentially expressed genes",
                      sets.bar.color = "#E69F00",
                      point.size = 5,
                      order.by = c("freq", "degree"),
                      show.numbers = "yes",
                      text.scale = c(2, 2, 1.3, 1.3, 2, 2))
  
  svglite::svglite(name_safe_svg, width = 14, height = 8)
  print(upset_plot)
  dev.off()
  
  cat("Saved upset as:", name_safe_svg, "\n")
}


so_upsetting <- map2(.x = biplot_list_in,
                     .y = de_direction, 
                     .f = f_upset_plot)
