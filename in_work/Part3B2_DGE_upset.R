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
library(UpSetR)

# 1.1 load def_variables.R and def_functions.R
source("/home/melanie/potato_RNAseq/RNAseq_potato_waterlogging/DGE_pipeline/def_variables.R")
source("/home/melanie/potato_RNAseq/RNAseq_potato_waterlogging/DGE_pipeline/def_functions.R")

# 1.2 Set wd
setwd(dir_base)

###################3
# generate tag for removing samples if defined
remove_samples <- f_gen_remove_samples(feature_columns)

# SET DIRECTION OF DEGs
DEG_DIRECTION <- c("DEGs", "up", "down")

# 1.4 in direcory
dir_dge <- paste0(
  dir_base, "DGE/"
)

set_name <- f_make_name(
  select_genes = select_genes,
  remove_samples = remove_samples
)

dir_in <- paste0(
  dir_base, 
  "DGE", "/",
  "consecutive_comps", "/",
  ref_gr, "_", 
  set_name
)

file_in_name <- paste(
  DEG_DIRECTION,
  "binema",
  paste0("padj-", pval),
  paste0("logFC-", log2FC),
  sep = "_"
)

dir_out <- dir_in

# 1.6 Check if in dir/file exist
if(dir.exists(dir_in)){
  cat()
}else{
  message("STOP! DEG directory does NOT! exist!. Check input directory")
}

if(file.exists(
  paste0(dir_in, "/",
         f_list_files(dir = dir_in,
                      match_str = file_in_name,
                      file_type = "csv"
                      )
)
)){
  cat("Filtered DGE Results exist. Continue creating binema...", "\n")
} else{
  message("STOP! Filterd DEG Results do NOT! exist! Check input requirements!")
}

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

names(biplot_list_in) <- paste0(DEG_DIRECTION , "_biplot")


f_upset_plot <- function(biplot_in, direction) {

name_safe_svg <- paste0(
  dir_out, "/", 
  paste(
    "Part3B2",
    direction,
    paste0("padj-", pval),
    paste0("logFC-", log2FC),
    "UPSET.svg",
    sep = "_"
    )
)

# Create the UpSet plot
upset_plot <- upset(data = biplot_in,
                  sets = rev(colnames(biplot_in)), # Customize based on desired order
                  nsets = length(colnames(biplot_in)),
                  keep.order = TRUE, # FALSE = orders sets by size
                  sets.x.label = "sum differentially expressed genes",
                  matrix.color = "black",
                  main.bar.color = "black",
                  mainbar.y.label = "number differentially expressed genes",
                  sets.bar.color = "black",
                  point.size = 5,
                  order.by = c("freq"), # options: c("freq", "degree")
                  show.numbers = "yes",
                  number.angles = 45,
                  text.scale = c(2, 2, 1.3, 1, 2, 2))

grid.text("My UpSet Plot Title", x = 0.5, y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))

svglite::svglite(name_safe_svg, width = 14, height = 8)
print(upset_plot)
dev.off()

cat("Saved upset as:", name_safe_svg, "\n")
}

##################################
f_upset_plot <- function(biplot_in, direction) {

name_safe_svg <- paste0(
  dir_out, "/", 
  paste(
    "Part3B2",
    direction,
    paste0("padj-", pval),
    paste0("logFC-", log2FC),
    "UPSET.svg",
    sep = "_"
  )
)

# Create the UpSet plot
upset_plot <- upset(data = biplot_in,
                  sets = rev(colnames(biplot_in)), # Customize based on desired order
                  nsets = length(colnames(biplot_in)),
                  keep.order = TRUE, # FALSE = orders sets by size
                  sets.x.label = "sum differentially expressed genes",
                  matrix.color = "black",
                  main.bar.color = "black",
                  mainbar.y.label = "number differentially expressed genes",
                  sets.bar.color = "black",
                  point.size = 5,
                  order.by = c("freq"), # options: c("freq", "degree")
                  show.numbers = "yes",
                  number.angles = 0,
                  text.scale = c(2, 2, 1.3, 0.1, 2, 2))

# Open the SVG device to save the plot
svglite::svglite(name_safe_svg, width = 20, height = 8)

# Print the plot to the SVG device
print(upset_plot)

# Add the title after printing the plot
grid::grid.text(
  paste( 
    direction,
    paste0("padj-", pval),
    paste0("logFC-", log2FC),
    sep = " "
  ),
  x = 0.5, 
  y = 0.98, 
  gp = grid::gpar(
    fontsize = 20, 
    fontface = "bold"
  )
)

# Close the SVG device
dev.off()

cat("Saved upset as:", name_safe_svg, "\n")
}



###################################
so_upsetting <- map2(.x = biplot_list_in,
                   .y = DEG_DIRECTION , 
                   .f = f_upset_plot)

