# Goal:   We might want to visualise filtered/unfiltered
#         expression of genes in a heatmap. This script generates a annotated,
#         unclustered and clustered heatmap, overlaying additional info on genes
# Input:  (filtered) list of DEGs, info for annotation of the heatmap 
# Output: heatmaps clustered, unclustered with annotation 

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

dir_in_mtx <- dir_goi

dir_out_de <- paste0(dir_goi, "/", name_out)

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

# 4. Load goi info
goi_file <- f_list_files(dir = paste0(dir_base, "DGE-goi"),
                         match_str = gois,
                         file_type = "csv")

goi_info <- f_read_csv(dir = paste0(dir_base, "DGE-goi", "/"),
                       csv_file = goi_file)

cat("number goi in genome:", nrow(goi_info))

# 5. Prep data for plotting annotated heatmap
# 5.1 log2 transform 
mtx_log <- log2(mtx+1)

# add info on chromosome
goi_info$chromosome <- sub(".*_(\\d{2})G.*", "\\1", rownames(goi_info))

# Get the indices where chromosome changes
gaps_row <- which(diff(as.numeric(goi_info$chromosome)) != 0)

# Remove goi not present in transcriptomics dataset
goi_info <- goi_info[rownames(goi_info) %in% rownames(mtx),]
cat("number goi in transcriptome:", nrow(goi_info))


# Ensure, all element in goi_info are factors
goi_info[] <- lapply(goi_info, as.factor)

# 5.2 Costume modify goi_info for annotated heatmap
#coldata_hm$`dn/ds` <- as.numeric(coldata_hm$`dn/ds`)

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

# "dn/ds" = c("#FFF5EE", colorRampPalette(c("#FFCF65", "#E69F00", "#8A5F00", "#362500"))(33)), #if continous scale

# 6.1 Draw a unclustered heat map
h2 <- pheatmap(mtx_log, 
               annotation_colors = my_colour,
               annotation_row = goi_info,
               cluster_cols = F,
               cutree_cols = F,
               cluster_rows = F,
               gaps_row = gaps_row
)

# 6.2 Safe unclustered heat map
ggsave(paste0(
  dir_out_de, "/",
  paste(
    "P3A2", heat_mtx, 
    "heatmap_unclustered_chromosomes", 
    "annotated",
    sep = "_", 
    ".svg"
  )
),
h2, 
device = "svg",
height = 8,
width =  13
)

dev.off()

# 7.1 Draw a clustered heat map
h3 <- pheatmap(mtx_log, 
               annotation_colors = my_colour,
               annotation_row = goi_info,
               cluster_cols = F,
               cutree_cols = F,
               cluster_rows = T
)

# 7.2 Safe clustered heat map
ggsave(paste0(
  dir_out_de, "/",
  paste(
    "P3A2", heat_mtx, 
    "heatmap_clustered", 
    "annotated",
    sep = "_", 
    ".svg"
  )
),
h3, 
device = "svg",
height = 8,
width =  13
)

dev.off()
