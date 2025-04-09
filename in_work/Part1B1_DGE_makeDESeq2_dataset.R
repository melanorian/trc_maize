# Goal: This script will create a DESeq2 dataset and DESeq normalised .csv count files from selected samples 
# Input: Requires a matrix input directory created in P1A2 with counts + metadata files
# Output: Directory with DESeq2 dataset with samples of interest + DESeqnormalised.csv countfiles

# 1.1 load libraries
library(DESeq2)
library(dplyr)

# 1.2 in/out directories

dir_in <- paste0(
  dir_base, "Input_Matrix/",
  f_make_name(
    select_genes = select_genes,
    remove_samples = remove_samples),
  "/"
  )

dir_out <- paste0(dir_base, "DDS/")
#dir_out_norm <- paste0(dir_base, "DESeq_norm/")

# ensure that dir_out exists

if (!file.exists(dir_out)) {
  dir.create(dir_out)
}

#if (!file.exists(dir_out_norm)) {
#  dir.create(dir_out_norm)
#}

# 2.Import files ----

# 2.1 Transcript level count matrix (CTS_trc) = cts, METADATA = coldata 
# Check if coldata and metadata are loaded and read files accordingly

f_load_if_absent_csv(
  dir_in = dir_in,
  name_env_variable = "coldata",
  match_str = "METADATA")

f_load_if_absent_csv(
  dir_in = dir_in,
  name_env_variable = "cts",
  match_str = "CTS")

# 3. Prepare for constructing DESeq2 object ----
# 3.1 round all counts to closest integer (DESeq doesn't accept decimals)

# 3.2 Define factor of interest
coldata <- f_define_factor(
  coldata = coldata,
  factor = factor_experiment
  )

# 3.3 Ensure factor_experiment is character
cat("Design formula:", as.character(design_formula), "\n")

# 4 Remove low read counts
cts_filtered <- eval(parse(text = filter_criterion))

# 4. Create dds datasets object ----
# Use the constructed formula in DESeqDataSetFromMatrix

dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData = coldata,
  design = as.formula(eval(parse(text = design_formula)))
  )

cat("constructed DESeq2 dataset (dds)...", "\n")

# 5. Extract normalised counts ----
# 5.1 Calculate size factors
dds <- estimateSizeFactors(dds)

# 5.2 Extract normalised counts
cts_norm <- counts(
  dds, 
  normalized = TRUE
  )

cat("generated DESeq2 normalised counts...")

# 6. Safe output ----

# 6.1 set name based on input
set_name <- f_make_name(
  select_genes = select_genes,
  remove_samples = remove_samples
  )

# 6.2 Safe output
# Create name for CTS_filtered based on butoff values
cutoff <- as.numeric(
  sub(".*cts >=([0-9]+).*", "\\1",
      filter_criterion)
  )
perc_samples<- as.numeric(
  sub(".*\\(([^/]+)/100\\).*", "\\1", 
      filter_criterion)
  )

name_cts_filtered <- paste(
  cutoff, 
  perc_samples, 
  sep = "_"
  )

write.csv(cts_filtered, 
          paste0(
            dir_in, "P1B1_" , 
          "cts_filtered", "_",
            name_cts_filtered,
            ".csv")
)

save(dds, 
     file = paste0(
       dir_out, "DDS_", 
       set_name, "_",
       "Des", design_shortname,
       ".RData")
     )

write.csv(cts_norm, 
          paste0(
            dir_in, "P1B1_" , 
            "DESeqNorm", ".csv")
          )

cat(
  "saved DDS file as:",
  paste0(
    dir_out, "DDS_", 
    set_name, "_",
    "Des", design_shortname,
    ".RData"),
  "\n"
  )

cat(
  "DESEq-normalised count matrix in:", 
  dir_in, "\n"
  )



rm(list = c("coldata", "cts", "cts_filtered", "cts_norm"))
