# Goal:   After the DGE analyis we might want to visualise the expression of a subest of DEG.
#         e.g. we are interested in the expression of DEG effectors. 
# Input:  List of DEGs, list of genes of interest (goi) (e.g. effecotr), count matrix or similar e.g.
#         DESeq2 normalised counts or CTS
# Output: A count matrix containing only the goi

# 1.1 clean environment and make sure all graphic devices are closed ----
rm(list = ls())

while (dev.cur() != 1) {
  dev.off()}

# 1.2 load libraries
library(dplyr)
library(tidyr)
library(readxl)
library(purrr)

# 1.2 Naming
# 1.1 load def_variables.R and def_functions.R
source("/media/wolf2/TKI2019/melanie/Code-transcriptomics-analysis/Pipeline/def_variables.R")
source("/media/wolf2/TKI2019/melanie/Code-transcriptomics-analysis/Pipeline/def_functions.R")

# 1.2 Set wd
setwd(dir_base)

# 1.3 in direcory

set_name <- f_make_name(select_genes = select_genes,
                         remove_samples = remove_samples)

# 1.4 out directory
dir_goi <- paste0(dir_base, "DGE-goi")

name_out <- set_name

dir_out <- paste0(dir_goi, "/", name_out, "_" ,gois)


if (!file.exists(dir_goi)) {
  dir.create(dir_goi)
}

if (!file.exists(dir_out)) {
  dir.create(dir_out)
}

# 2. Import Matrixes ----
# 2.1 If "MTX-name"  in plot_mtx load it

if ("CTS" %in% plot_mtx) {
  # Construct file path to cts
  dir_in_cts <- paste0(dir_base, "Input_Matrix/",
                       set_name,
                       "/")
  # If "CTS" is present, then execute the read function
  cts_in <- f_list_files(dir = dir_in_cts,
                         match_str = "CTS",
                         file_type = "csv")
  
  cts <- f_read_csv(dir = dir_in_cts,
             csv_file = cts_in)
}


if ("CPM" %in% plot_mtx) {
  # Construct file path to cts
  dir_in_cpm <- paste0(dir_base, "Input_Matrix/",
                       set_name,
                       "/")
  # If "CTS" is present, then execute the read function
  cpm_in <- f_list_files(dir = dir_in_cpm,
                         match_str = "CPM",
                         file_type = "csv")
  
  cpm <- f_read_csv(dir = dir_in_cpm,
                    csv_file = cpm_in)
}

if ("DESeq" %in% plot_mtx) {
  # Construct file path to cts
  dir_in_des <- paste0(dir_base, "Input_Matrix/",
                       set_name,
                       "/")
  # If "CTS" is present, then execute the read function
  des_in <- f_list_files(dir = dir_in_des,
                         match_str = "DESeqNorm",
                         file_type = "csv")
  
  des <- f_read_csv(dir = dir_in_des,
                    csv_file = des_in)
}

# 2.2 If "MTX-name" in plot_mtx add to list_mtx and remove individual mtx
# Get all objects in the global environment
objects <- ls()

# Create an empty list to store the lists
list_mtx <- list()

# Check if 'cts' is present and is a list
if ("cts" %in% objects && is.list(get("cts"))) {
  list_mtx <- c(list_mtx, list(cts = get("cts")))
  rm(cts)
}

# Check if 'cpm' is present and is a list
if ("cpm" %in% objects && is.list(get("cpm"))) {
  list_mtx <- c(list_mtx, list(cpm = get("cpm")))
  rm(cpm)
}

# Check if 'des' is present and is a list
if ("des" %in% objects && is.list(get("des"))) {
  list_mtx <- c(list_mtx, list(des = get("des")))
  rm(des)
}

names_mtx <- names(list_mtx)
cat("Created Vector:", names_mtx)

# 3. Import list genes of interest
# 3.6 load list of DEGs
goi <- read_excel(goi_path)
goi <- goi$Column1
goi <- sub("-T1$", "", goi)

cat("loaded goi file. Name samples:", head(goi))
cat("nr genes in list:", length(goi))

# 5. Filter by goi
# function defined in def_functions.R
list_mtx_filtered <- lapply(list_mtx, f_filter_df, goi_names = goi)

# Print the structure of the filtered list
cat("Names of list_mtx_filtered:")
head(names(list_mtx_filtered))

# 6. Safe resuls 
# Apply the function to save each dataframe in list_mtx_filtered
map2(list_mtx_filtered, names_mtx, ~ f_safe_df(df = .x, name_df = .y, dir = dir_out, gois = gois, part = "P1A4"))

cat("saved count matrix filtered for goi at:", dir_out)
