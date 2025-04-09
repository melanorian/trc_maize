# Goal:   This script allows creating a binary matrix for DGE-GOI data in order to 
#         plot the results as an upsetplot, and identify which genes are uniquely expressed
#         at a specific timepoint or over a set number of timepoints
# Input:  results of DEG analysis in .RData format where res$contrast$direction 
# Output: Binary matrix, upsetplot, list of genes uniquely DE at a specific timepoint or 
#         expressed over all timepoints (possibly interesting candidates)

# extract differentially expressed genes
# clean environment and make sure all graphic devices are closed

# 1. Set-up & import data ----
# Clean environment
rm(list = ls())

while (dev.cur() != 1) {
  dev.off()}

# 1.1 load libraries
library(dplyr)
library(purrr)

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
# 2.1 Construct file name based on input
file_in_name <- paste("filterLFC", paste0("FC", log2FC), 
                      paste0("PADJ", pval), sep = "_"
)

# 2.2 Identify file matching file_in
file_in <- f_list_files(dir = dir_in,
                        match_str = "DEG-filtered",
                        file_type = "RData")

# 2.3 construct path to file
path_in <- paste0(dir_in, "/", file_in)

# 2.4 load file
DEG <- f_assign_load(dir = path_in)

# 3. Contrasts analysed in DEG
list_contr <- names(DEG)

# 4. Seperate sets of DEG (all DEGs, up or down regulated)
# 4.1 Function to extract goi from all contrasts
f_select_DEGs <- function(contrast, direction){
  gene_df <- DEG[[contrast]][[direction]]
  gene_list <- rownames(gene_df)
  return(gene_list)
}

# 4.2 extract list of all DEG
for (direction in de_direction){
  # dynamically create a list name based on available DE_directions
  variable_name <- paste0(direction, "_list")
  assign(variable_name, list())
  # Generate the lists of gene names differentially  expressed at each timepoint
  result <- lapply(list_contr, f_select_DEGs, direction = direction)
  # Ensure the naming of list elements is maintained
  names(result) <- list_contr
  # store the results in the list with the assigned name
  assign(variable_name, result)
}

# 4.3 Make a list of name lists so that next functions can be run over all
de_list <- c(mget(paste0(de_direction, "_list")))

# 5. Create a matrix with rows corresponding to DE of any contrast
# columns to contrasts. 

# 5.1 For loop to crate the mtx for each direction of comparison
for(direction in de_direction){
  # Dynamically name empty matrix based on available DE_directions
  mtx_name <- paste0(direction, "_mtx")
  # Identify all DEG gene names for assigning rownames in mtx
  unlist <- de_list[[paste0(direction, "_list")]] %>% unlist() %>% unique()
  
  # Initialize an empty matrix with proper row and column names
  mtx_result <- matrix(NA, nrow = length(unlist), 
                       ncol = length(list_contr), 
                       dimnames = list(unlist, list_contr))
  
  # Fill the matrix with gene names when they are DE in a given contrast
  for (contr in 1:length(list_contr)) {
    present_genes <- de_list[[paste0(direction, "_list")]][[contr]]
    mtx_result[unlist %in% present_genes, contr] <- unlist[unlist %in% present_genes]
  }
  
  # Convert the matrix to a data frame
  mtx_result <- as.data.frame(mtx_result)
  
  # Assign mtx_result to the dynamically created matrix
  assign(mtx_name, mtx_result)
}

cat("generated mtx...")


# 5.2 Make a list of mtx so that next functions can be run over all sets
mtx_list <- c(mget(paste0(de_direction, "_mtx")))

# 6. For some purposes e.g. upset plot convert into binary matrix
# 6.1 For loop to convert each matrix (dim: row = DEGs, col = contr)
for (direction in de_direction){
  # Dynamically name empty matrix based on available DE_directions
  binema_name <- paste0(direction, "_binema")
  
  binary_df <- apply(mtx_list[[paste0(direction, "_mtx")]], 2, 
                     function(x) ifelse(is.na(x), 0, 1))
  binary_df <- as.data.frame(binary_df)
  assign(binema_name, binary_df)
  }

# 6.2 store res in vector
mtx_binema <- c(mget(paste0(de_direction, "_binema")))


# 7. Save output
# 7.1 Matrix character 
map2(.x = mtx_binema,
     .y = paste0(de_direction, "_mtx"),
     .f = f_safe_df,
     dir = dir_out, 
     gois = "DE-character",
     part = "P3B")

# 7.2 Binary matrix
map2(.x = mtx_binema,
     .y = paste0(de_direction, "_binema"),
     .f = f_safe_df,
     dir = dir_out, 
     gois = "",
     part = "P3B")
