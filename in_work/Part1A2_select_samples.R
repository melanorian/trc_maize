# Goal: This script allows to drop genesets from the RNAseq dataset
# Input: Requires a matrix with counts + metadata file 
# Output: Directory with counttable.csv, METADATA files and samples of interest

dir_in_1A2 <- paste0(
  dir_base, 
  "Input_Matrix/",
  f_make_name(select_genes = select_genes,
              remove_samples = remove_samples
              )
)

if(dir.exists(dir_in_1A2)){
  cat("input for sample selection exists, continue..")
} else {
  cat("Input directory with removed samples does NOT! exist")
}

dir_out <- paste0(dir_base, "Input_Matrix/",
                  f_make_name(select_genes = select_genes,
                              remove_samples = remove_samples),
                  "/"
                  )



# create output directory if it doesn't exist
if (!file.exists(dir_out)) {
  dir.create(dir_out)
}

# 2.Import files ----

# 2.1 Transcript level count matrix (CTS_trc) = cts 
cts_file <- f_list_files(dir = dir_in_1A2, match_str = "CTS", file_type = "csv")
cts <- f_read_csv(dir = dir_in_1A2, csv_file = cts_file)

cat("loaded cts...", "\n")

# 2.2 METADATA = coldata
meta_file <- f_list_files(dir = dir_in_1A2, match_str = "METADATA", file_type = "csv")
coldata <- f_read_csv(dir = dir_in_1A2, csv_file = meta_file)
cat("loaded metadata...", "\n")

# 3. Filter for selected genes defined in def_variables.R

# Check if select_genes is NULL
if (is.null(select_genes)) {
  cts_keep <- cts # If select_genes is NULL, keep all rows
} else {
  # If select_genes is not NULL, subset rows based on the pattern
  cts_keep <- cts[grepl(select_genes, rownames(cts)), , drop = FALSE]
}

cat("Select for further analysis", nrow(cts_keep), " of ", nrow(cts), "genes present in dataset", "\n")
rm(cts)

# 4. Save output
write.csv(cts_keep, paste0(dir_out, "/", "P1A2_CTS_trc", ".csv"))
write.csv(coldata, paste0(dir_out, "/", "P1A2_METADATA", ".csv"))

cat("finished saving files in:", dir_out, "\n")
