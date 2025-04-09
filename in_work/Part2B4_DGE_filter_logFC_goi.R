# Goal:   Filter matrices of DEG results of interest e.g. log2FC for goi
# Input:  Matrix of with DEG results e.g. LFC
# Output: A filtered DEG list and a list of gene names in the filtered list

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

# 1.3 Load list goi
goi <- read_excel(goi_path)
goi <- goi$Column1
goi <- sub("-T1$", "", goi)

cat("loaded goi file. Name samples:", head(goi))
cat("nr genes in list:", length(goi))

# 1.4 in direcory

set_name <- f_make_name(select_genes = select_genes,
                        remove_samples = remove_samples)

dir_in <- paste0(dir_base, "DGE", "/",
                 ref_gr, "_", set_name
)

# 1.5 out directory
name_out <- paste(paste0("FC", log2FC), 
                  paste0("PADJ", pval), 
                  sep = "_"
)

dir_de_goi <- paste0(dir_base, "DGE-goi", "/",
                 set_name, "_", gois
)

de_parameters <- paste("unfilter", paste0("FC", log2FC), 
                      paste0("PADJ", pval), sep = "_"
)

dir_out <- paste0(dir_de_goi, "/", ref_gr, "_", de_parameters)

if (!file.exists(dir_out)) {
  dir.create(dir_out)
}

# 2.Import mtx with DEG results
# 2.1 Construct file name based on input
file_in_LFC <- f_list_files(dir = dir_in,
                        match_str = "LFC-ALL",
                        file_type = "csv")

file_in_padj <- f_list_files(dir = dir_in,
                            match_str = "padj-ALL",
                            file_type = "csv")

# 2.3 Identify file matching file_in

# 2.4 construct path to file
path_in_LFC <- paste0(dir_in, "/", file_in_LFC)
path_in_padj <- paste0(dir_in, "/", file_in_padj)

# 2.5 load file
lfc <- f_read_csv(dir = paste0(dir_in, "/"), file_in_LFC)
padj <- f_read_csv(dir = paste0(dir_in, "/"), file_in_padj)

# 3. Filter mtx for effectors
logFC_effectors <- lfc[row.names(lfc) %in% goi, , drop = FALSE]
padj_effectors <- padj[row.names(padj) %in% goi, , drop = FALSE ]

# 4. Safe results 
write.csv(logFC_effectors, file = paste0(dir_out, "/", "P2B_", "LFC-", gois, ".csv"), row.names = TRUE)
write.csv(padj_effectors, file = paste0(dir_out, "/", "P2B_", "padj-", gois, ".csv"), row.names = TRUE)

cat("saved results in: ", dir_out)
