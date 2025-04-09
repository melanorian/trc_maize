# Goal: This plot generates a volcano plot to visualise DGE
# Input: .csv files of log2FC from DGE analysis and adjusted p-values
# Output: volcano plot as .png (.svg is larger, change file type if needed)

# clean environment and make sure all graphic devices are closed
rm(list = ls())

while (dev.cur() != 1) {
  dev.off()}

# 1. base set up install/load libraries, setwd, read files ----
# 1.1 load libraries
#BiocManager::install('EnhancedVolcano')
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)
library(svglite)
library(gridExtra)

# 1.2 Naming
# 1.1 load def_variables.R and def_functions.R
source("/home/melanie/potato_RNAseq/RNAseq_potato_waterlogging/DGE_pipeline/def_variables.R")
source("/home/melanie/potato_RNAseq/RNAseq_potato_waterlogging/DGE_pipeline/def_functions.R")

# 1.2 Set wd
setwd(dir_base)

# generate tag for removing samples if defined
remove_samples <- f_gen_remove_samples(feature_columns)

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

# file_in <- paste0(
#     paste(
#     "P2A", "unfiltered", 
#     paste0("FC", "default"), 
#     paste0("PADJ", pval),
#     paste0("Des", design_shortname),
#     sep = "_"),
#     ".RData"
#   )

file_in <- "DGE_RESULTS-joined.RData"

# 1.5 out directory
dir_out <- dir_in

name_out <- paste(
  paste0("FC", log2FC),
  paste0("PADJ", pval), 
  sep = "_"
  )

# 1.6 Check if in dir/file exist
if(dir.exists(dir_in)){
  cat()
}else{
  message("STOP! DEG directory does NOT! exist!. Check input directory")
}

if(file.exists(
  paste0(
    dir_in, "/",
    file_in)
  )){
  cat("DGE Results exist. Continue volcano plotting...", "\n")
} else{
  message("STOP! DEG Results do NOT! exist! Check input requirements!")
}

# 2. Import DEG input (file with all results unfiltered)
# 2.1 Load file with all results of DGE
DEG <- f_assign_load(
  dir = paste0(
    dir_in, "/",
    file_in
  )
  )

# 2.4 remove NA if needed
DEG <- lapply(DEG, na.omit)
cat("omitted NAs")

# 3. Set up for plotting
# 3.1 what should be plotted on the Axis?
x_axis <- 'log2FoldChange'
y_axis <- 'padj'

cat("Plotted on volcano x-axis", x_axis)
cat("Plotted on volcano y-axis", y_axis)

# 3.2 Dimensions of grid
grid_rows <- round(sqrt(length(DEG))) # how many rows should the grid have
grid_cols <- round(sqrt(length(DEG))) # how many columns should the grid have

# 3.3 Extract list contrasts for accessing data 
list_contr <- names(DEG)
cat("List for plotting volcano:", list_contr)

# 4. Main function for Volcano plot ----
# https://github.com/kevinblighe/EnhancedVolcano

#    Functions to draw volcano plots from a list,
#    create a common figure legend
#    arrange everything in a grid panel

f_explode_volcano <- function(DEG_list, x_axis, y_axis, pval, log2FC, grid_rows, grid_cols) {
  # Function to draw the volcano plot
  f_volcano <- function(DEG, x_axis, y_axis, pval, log2FC, title) {
  volcano_plot <- EnhancedVolcano(DEG,
                                  lab = rownames(DEG),
                                  x = x_axis,
                                  y = y_axis,
                                  pCutoff = pval,
                                  FCcutoff = log2FC,
                                  pCutoffCol = "padj",
                                  cutoffLineType = "twodash",
                                  cutoffLineCol = "grey",
                                  cutoffLineWidth = 0.8,
                                  max.overlaps = 15,
                                  gridlines.major = FALSE,
                                  gridlines.minor = FALSE,
                                  ylim = c(0, 300),
                                  pointSize = 2,
                                  labSize = 2.0,
                                  col=c("#999999", "#E69F00", "#56B4E9", "#009E73"),
                                  title = title,
                                  titleLabSize = 10,
                                  subtitle = NULL,
                                  axisLabSize = 10,
                                  #drawConnectors = TRUE,
                                  widthConnectors = 0.1)
    return(volcano_plot)
  }
  
  # Function to extract figure legend from plots
  f_get_legend <- function(myggplot) {
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }
  
  # Function to draw a grid panel of volcano plots
  f_draw_grid <- function(volcano_list, grid_rows, grid_cols, legends) {
    gridExtra::grid.arrange(grobs = lapply(seq_along(volcano_list), function(i) {
      p <- volcano_list[[i]]
      legend <- legends[[i]]
      p + theme(legend.position = "none")  # Remove individual legends
    }), nrow = grid_rows, ncol = grid_cols)
  }
  
  # Draw volcano plots
  volcano_list <- lapply(seq_along(DEG_list), function(i) {
    DEG <- DEG_list[[i]]
    title <- names(DEG_list)[i]
    f_volcano(DEG, x_axis, y_axis, pval, log2FC, title)
  })
  
  # Extract legends
  legends <- lapply(volcano_list, f_get_legend)
  
  # Draw grid panel
  grid_volcano <- f_draw_grid(volcano_list, grid_rows, grid_cols, legends)
  
  return(grid_volcano)
}

# 5. Draw plots
# 5.1 Draw Plot for goi
volcano <- f_explode_volcano(DEG_list = DEG,
                                 x_axis = x_axis,
                                 y_axis = y_axis,
                                 pval = pval,
                                 log2FC = log2FC,
                                 grid_rows = grid_rows,
                                 grid_cols = grid_cols)

# 6.1 Save plot all genes
ggsave(paste0(
  dir_out, "/",
  paste(name_out, "volcano_all", sep = "_"), ".png"),
  volcano,
  device = "png",
  width = 15,
  height = 15
  )

cat("volcano plot for all genes saved at:", dir_out)

# # Safe individual plots all
# f_safe_one_volcano <- function(volcano_list, contr, dir_out, name_plot){
#   one_volcano <- volcano_list[[contr]]
#   print(paste0("printing:", contr))
#   ggsave(paste0(dir_out, "/", paste(contr, name_plot, "labeled", sep = "_"), ".png"), 
#          one_volcano, 
#          device = "png",
#          width = 25,
#          height = 20)
#   print(paste0("saved:", contr))
# }
# 
# 
# lapply(list_contr, f_safe_one_volcano, 
#        volcano_list = volcano,
#        dir_out = dir_out_all,
#        name_plot = "all")