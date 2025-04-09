# Goal:   After the DGE analyis we might want to visualise filtered/unfiltered
#         results of DGE analysis in a heatmap. This script generates a simple,
#         unclustered heatmap
# Input:  (filtered) expression matrix
# Output: heatmap

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

set_name <- f_make_name(select_genes = select_genes,
                        remove_samples = remove_samples)

# 1.5 out directory
dir_goi <- paste0(dir_base, "DGE-goi", "/", 
                  set_name, "_", gois
)

name_out <- paste(ref_gr, "filterLFC", 
                  paste0("FC", log2FC), 
                  paste0("PADJ", pval), 
                  sep = "_"
)

dir_in_de <- paste0(dir_goi, "/", name_out)
dir_in_mtx <- dir_goi

# 1.6 Define out directory
dir_out <- dir_in_de

# 2. Load DEG files
# 2.1 Identify file matching file_in
file_in_de <- f_list_files(dir = dir_in_de,
                        match_str = file_in_name,
                        file_type = "RData")

# 2.2 construct path to file
path_in <- paste0(dir_in_de, "/", file_in_de)

# 2.3 load file
DEG <- f_assign_load(dir = path_in)

# manual import for all genes
DEG <- f_assign_load(dir = "/media/wolf2/TKI2019/melanie/Salmon_decoy_ISR_utr500_CDS/DGE/M_genes-S_samples-no-S-GS/P2A_filterLFC_FC1.5_PADJ0.05.RData")

# 3. Load mtx
file_in_de <- f_list_files(dir = dir_in_mtx,
                           match_str = heat_mtx,
                           file_type = "csv")


if (length(file_in_de) == 1) {
  dir_in_mtx <- paste(dir_in_mtx, file_in_de, sep = "/")
  # If "CTS" is present, then execute the read function
  mtx_in <- f_list_files(dir = dir_in_mtx, match_str = heat_mtx, file_type = "csv")
  mtx <- f_read_csv(dir = dir_in_mtx, csv_file = mtx_in)
  cat("count file loaded...")
} else {
  message("no or >1 mtx matching the conditions!")
}

# manual import cts for all genes
# mtx <- f_read_csv(dir = "/media/wolf2/TKI2019/melanie/Salmon_decoy_ISR_utr500_CDS/Input_Matrix/genes-S_samples-no-S-GS/",
#             csv_file = "P1A2_CTS_trc.csv")

# 4. Load goi info
# 3.6 load list of DEGs
goi <- read_excel(goi_path)
names(goi)[1] <- "goi_ID"
goi$goi_ID <- sub("-T1$", "", goi$goi_ID)

# 5. Plot heatmap
# 5.1 log2 transform 
mtx_log <- log2(mtx+1)

# 5.2 define colours
paletteLength <- 50
myColor <- colorRampPalette(c("#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF", "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026"))(paletteLength)


h1 <- pheatmap(mtx_log,
               color = myColor,
               cluster_cols = FALSE,
               cluster_rows = F)

ggsave(paste0(
  dir_out, "/",
  paste(
    "P3A1", heat_mtx, 
    "heatmap_clustered", 
    sep = "_", 
    ".svg"
    )
  ),
  h1, 
  device = "svg",
  height = 8,
  width =  10
  )

dev.off()

cat("Heatmap of gene expression saved at:", dir_out)

h2 <- pheatmap(mtx_log,
               color = myColor,
               cluster_cols = FALSE,
               cluster_rows = T)

ggsave(paste0(
  dir_out, "/",
  paste(
    "P3A1", heat_mtx, 
    "heatmap_clustered", 
    sep = "_", 
    ".svg"
  )
),
h2, 
device = "svg",
height = 8,
width =  10
)

dev.off()

cat("Heatmap of gene expression saved at:", dir_out)

