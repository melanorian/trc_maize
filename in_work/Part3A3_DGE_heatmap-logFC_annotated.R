# Goal:   Visualise logFC2 as heatmap
#         unclustered heatmap, overlaying additonal info on genes
# Input:  here we use the full list of DEGs, possibly info with additional annotation 
# Output: heatmaps with log2FC

# 1.1 clean environment and make sure all graphic devices are closed ----
rm(list = ls())

while (dev.cur() != 1) {
  dev.off()}

# 1.2 load libraries
library(readxl)

# 1.2 Naming
# 1.1 load def_variables.R and def_functions.R
source("/media/wolf2/TKI2019/melanie/Code-transcriptomics-analysis/Pipeline/def_variables.R")
source("/media/wolf2/TKI2019/melanie/Code-transcriptomics-analysis/Pipeline/def_functions.R")

# 1.2 Set wd
setwd(dir_base)

# 1.4 in direcory
base_in <- 

set_name <- f_make_name(select_genes = select_genes,
                        remove_samples = remove_samples)

# 1.5 out directory
dir_DGE_goi <- paste0(dir_base, "DGE-goi", "/", 
                  set_name, "_", gois
)

dir_out_name <- paste(ref_gr, "unfilter", 
                     paste0("FC", log2FC), 
                     paste0("PADJ", pval), 
                     sep = "_"
)

dir_out <- paste0(dir_DGE_goi, "/", dir_out_name)

dir_in <- dir_out

# 2. Load DEG files
# 2.1 Identify file matching file_in
file_in_lfc <- f_list_files(dir = dir_in,
                           match_str = "LFC",
                           file_type = "csv")

file_in_padj <- f_list_files(dir = dir_in,
                            match_str = "padj",
                            file_type = "csv")

# 2.3 load file
lfc <- f_read_csv(dir = dir_in, paste0("/", file_in_lfc)) 
padj <- f_read_csv(dir = dir_in, paste0("/", file_in_padj)) 

# 4. Load goi info
goi_file <- f_list_files(dir = paste0(dir_base, "DGE-goi"),
                         match_str = gois,
                         file_type = "csv")

goi_info <- f_read_csv(dir = paste0(dir_base, "DGE-goi", "/"),
                       csv_file = goi_file)

cat("number goi in genome:", nrow(goi_info))

# Remove goi not present in transcriptomics dataset
goi_info <- goi_info[rownames(goi_info) %in% rownames(lfc),]
cat("number goi in transcriptome:", nrow(goi_info))

# Ensure, all element in goi_info are factors
goi_info[] <- lapply(goi_info, as.factor)

# 5. Prep data for plotting annotated heatmap
# 5.3 Define costum colours for annotation based on goi_info
my_colour <- list(
  variation = c(core = "#999999", accessory = "#E69F00", unique = "#D55E00"),
  type = c(na = "#FFF5EE", RXLR = "#56B4E9", CRN = "#F0E442"),
  selection = c(na = "#FFF5EE", positive = "#E69F00", stabilising = "#009E73"),
  cluster = c(clustered = "#0072B2", unclustered = "#D55E00"),
  motif = c("na" = "#FFF5EE", "EER" = "#88CCEE", "EER,RXLR" = "#CC6677",
            "LFLAK" = "#DDCC77", "LFLAK,DWL" = "#117733", "LFLAK,DWL,HVLV" = "#332288",
            "HVLV" = "#AA4499", "RXLR" = "#44AA99", "RXLR,EER" = "#999933",
            "WY,RXLR" = "#882255", "WY,RXLR,EER" = "#661100"
  )
)

# 3.2 Save filtered cts file
#save_mtx <- padj
#name_output_logFC <- paste(pre, part, "PADJ", sep = "_")
#write.csv(save_mtx, paste0(name_output_logFC, "_effectors.csv"), row.names = TRUE)

# transform mtx with padj into symbols to indicate significance
thresholds <- c(0.05, 0.01, 0.001) # Define the thresholds and corresponding symbols
symbols <- c("x", "xx", "xxx")

mtx_sym_padj <- matrix("", nrow = nrow(padj), ncol = ncol(padj),
                       dimnames = list(rownames(padj), colnames(padj)))

# Replace values in the new matrix according to thresholds
for (i in 1:length(thresholds)) {
  mtx_sym_padj[padj <= thresholds[i]] <- symbols[i]
}

# Draw heatmap with log2FC results
min_val <- min(lfc[!is.na(lfc)])
max_val <- max(lfc[!is.na(lfc)])

# Define the number of steps below -1.5
breaks_sign_neg <- c(round(min_val):-1.5)
colour_neg <- colorRampPalette(c("#313695", "#4575B4", "#74ADD1", "#ABD9E9"))(length(breaks_sign_neg))

# Define additional breaks around zero for bright colors
breaks_ns <- seq(-1.4, 1.4, length.out = 5)
colour_ns <- colorRampPalette(c("#E0F3F8", "white", "#FFFFBF"))(length(breaks_ns))

# Define the number of steps above -1.5
breaks_sign_pos <- c(1.5:round(max_val))
colour_pos <- colorRampPalette(c("#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026"))(length(breaks_sign_pos))

# Summarize for pheatmap
myBreaks <- c(breaks_sign_neg, breaks_ns, breaks_sign_pos)
myColor <- c(colour_neg, colour_ns, colour_pos)

# does not accept 0s
lfc[is.na(lfc)] <- 0

# Draw Unclustered heatmap
h3 <- pheatmap(lfc,
               color = myColor,
               cluster_cols = FALSE,
               cluster_rows = FALSE,
               breaks = myBreaks,
               na_col = "grey",
               annotation_colors = my_colour,
               annotation_row = goi_info,
               display_numbers = mtx_sym_padj,
               number_color = "black",
               angle_col = 45)

ggsave(paste0(
  dir_out, "/",
  paste(
    "P3A3", heat_mtx, 
    "heatmap_lfc", 
    "annotated_unclustered",
    sep = "_", 
    ".svg"
  )
),
h3, 
device = "svg",
height = 45,
width =  13
)

dev.off()

# Draw clustered heatmap
h4 <- pheatmap(lfc,
               color = myColor,
               cluster_cols = FALSE,
               cluster_rows = TRUE,
               breaks = myBreaks,
               na_col = "grey",
               annotation_colors = my_colour,
               annotation_row = goi_info,
               display_numbers = mtx_sym_padj,
               number_color = "black",
               angle_col = 45)

ggsave(paste0(
  dir_out, "/",
  paste(
    "P3A3", heat_mtx, 
    "heatmap_lfc", 
    "annotated_clustered",
    sep = "_", 
    ".svg"
  )
),
h4, 
device = "svg",
height = 45,
width =  13
)

dev.off()

cat("results saved at:", dir_out)
