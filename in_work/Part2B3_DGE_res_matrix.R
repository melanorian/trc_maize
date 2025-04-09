# Goal:   After the DGE analyis we might want to visualise/analyse results by 
#         analysing log2FC/and or padj. Here we can make matrixed of log2FC
# Input:  Unfiltered output of DGE, list of genes of interest to filter against
#         usually indicated as "filterLFC" in /DEG directory
# Output: A filtered DEG list and a list of gene names in the filtered list

# 1.1 clean environment and make sure all graphic devices are closed ----
rm(list = ls())

while (dev.cur() != 1) {
  dev.off()}

# 1.2 load libraries
library(dplyr)
library(tidyr)
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

# 1.5 out directory
dir_in <- paste0(dir_base, "DGE", "/",
                 ref_gr, "_", set_name
)

name_out <- paste(paste0("FC", log2FC), 
                  paste0("PADJ", pval), 
                  sep = "_"
)

dir_out <- dir_in

# 2.1 Import DEG results


# 2.2 Construct file name based on input
file_in_name <- paste("filterLFC", paste0("FC", log2FC), 
                      paste0("PADJ", pval), sep = "_"
)

# 2.3 Identify file matching file_in
file_in <- f_list_files(dir = dir_in,
                        match_str = "unfiltered",
                        file_type = "RData")

# 2.4 construct path to file
path_in <- paste0(dir_in, "/", file_in)

# 2.5 load file
DEG <- f_assign_load(dir = path_in)

# 3. filter DEG list for genes of interest 
# 3.1 list the contrasts in DEG
list_contrasts <- names(DEG)
list_contrasts <- list_contrasts[list_contrasts != "Intercept"]

# 3.2 use f_filter_df from def_functions.R to filter each element in the lists
DEG_filtered <- lapply(list_contrasts, function(contrast) {
  # Get the DEG object for this contrast
  contrast_df <- DEG[[contrast]]
  df_filtered <- f_filter_df(contrast_df, goi_names = goi) # def in def_functions.R
  return(df_filtered)
})

# 3.3 Name filtered list with contrast names
names(DEG_filtered) <- list_contrasts

# 4. Filter results for value of interest, usually LFC
# 4.1 Define a function to return a matrix of values of interest from DEGs
f_DEG_mtx <- function(DEG, contr, mtx_values){ # can be: "log2FoldChange", padj,..
  
  # create an empty matrix rows = genenames, columns = contrast
  genes_name <- rownames(DEG[[contr[1]]]) # uses first element in contrast as reference
  print(genes_name)
  mtx <- matrix(NA, nrow = length(genes_name), ncol = length(contr),
                dimnames = list(genes_name, contr))
  
  # Iterate over each contrast
  for (contrast_index in seq_along(contr)) { # use seq_along to access the indices of contrast instead of elements 
    current_contrast <- contr[contrast_index]
    # Extract the specified matrix values for the current contrast
    values <- DEG[[current_contrast]][[mtx_values]]
    print(values)
    # Fill the corresponding column in the matrix
    mtx[, contrast_index] <- values
  }
  
  return(mtx)
}

# 4.2 Use function to generate matrixes with log2FC and adjusted pvalues
mtx_log <- f_DEG_mtx(DEG = DEG, contr = list_contrasts, mtx_values = "log2FoldChange")
mtx_padj <- f_DEG_mtx(DEG = DEG, contr = list_contrasts, mtx_values = "padj")

# 5. Safe results 
write.csv(mtx_log, file = paste0(dir_out, "/", "P2B_", "LFC-ALL", ".csv"), row.names = TRUE)
write.csv(mtx_padj, file = paste0(dir_out, "/", "P2B_", "padj-ALL", ".csv"), row.names = TRUE)

cat("saved results in: ", dir_out)
