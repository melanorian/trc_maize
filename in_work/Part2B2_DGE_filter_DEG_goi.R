# Goal:   After the DGE analyis we might want to filter results for specific
#         genes of interest e.g. effectors for downstream analysis
# Input:  List of (significant) DEGs, list of genes of interest (goi)
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
dir_goi <- paste0(dir_base, "DGE-goi", "/", 
                  set_name, "_", gois
                  )

name_out <- paste(ref_gr, "filterLFC", 
                  paste0("FC", log2FC), 
                  paste0("PADJ", pval), 
                  sep = "_"
                  )

dir_out <- paste0(dir_goi, "/", name_out)


if (!file.exists(dir_goi)) {
  dir.create(dir_goi)
}

if (!file.exists(dir_out)) {
  dir.create(dir_out)
}

# 2.1 Import DEG results
dir_in <- paste0(dir_base, "DGE", "/",
                 ref_gr, "_", set_name
                 )

# 2.2 Construct file name based on input
file_in_name <- paste("filterLFC", paste0("FC", log2FC), 
                 paste0("PADJ", pval), sep = "_"
                 )

# 2.3 Identify file matching file_in
file_in <- f_list_files(dir = dir_in,
                        match_str = file_in_name,
                        file_type = "RData")

# 2.4 construct path to file
path_in <- paste0(dir_in, "/", file_in)

# 2.5 load file
DEG <- f_assign_load(dir = path_in)

# 3. filter DEG list for genes of interest 
# 3.1 list the contrasts in DEG
list_contrasts <- names(DEG)

# 3.2 use f_filter_df from def_functions.R to filter each element in the lists
DEG_filtered <- lapply(list_contrasts, function(contrast_name) {
  # Get the DEG object for this contrast
  contrast <- DEG[[contrast_name]]
  # Iterate over each direction in the contrast
  filtered_directions <- lapply(contrast, function(direction) {
    f_filter_df(direction, goi_names = goi)
  })
  names(filtered_directions) <- names(contrast)
  return(filtered_directions)
})

# 3.3 Name outer list with contrast names
names(DEG_filtered) <- list_contrasts

# 4. Extract gene names present in each list
# 4.1 Function to extract gene names from all contrasts
f_deg_names <- function(contrast, direction, DEG_list){
  gene_df <- DEG_list[[contrast]][[direction]]
  gene_list <- rownames(gene_df)
  return(gene_list)
}

# 4.2 Define function to extract gene names from all contrasts in a DEG set
get_gene_names <- function(de_direction, list_contrasts, DEG_filtered) {
  # Apply the function for each direction
  names_DEG_list <- lapply(de_direction, function(direction) {
    # Apply the function to each contrast and store the result in a list
    result <- lapply(list_contrasts, 
                     f_deg_names, 
                     direction = direction,
                     DEG_list = DEG_filtered)
    # Name the resulting list with contrast names
    names(result) <- list_contrasts
    return(result)
  })
  return(names_DEG_list)
}

# 4.3 Call the function
names_DEG_list <- get_gene_names(de_direction, list_contrasts, DEG_filtered)
names(names_DEG_list) <- de_direction 

# 4.4 Print the structure of names_DEG_list
print(str(names_DEG_list))

# 5. Safe output
# 5.1 safe filtered list of DEG
save(DEG_filtered, file = paste0(dir_out, "/", "P2B_DEG-filtered", ".RData"))
save(names_DEG_list, file = paste0(dir_out, "/", "P2B_DEG-names-filtered", ".RData"))

cat("saved output at:", dir_out)
