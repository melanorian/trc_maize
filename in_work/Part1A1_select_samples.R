# Goal: This script allows to drop samples from the RNAseq dataset
# Input: Requires a matrix with counts + metadata file 
# Output: matching METADATA and CTS (count matrix).csv files with the samples of interest

dir_in <- paste0(dir_base, "Input_Matrix/genes-all_samples-all/") # genes-all_samples-all generated in P0

if(dir.exists(dir_in)){
  cat("Input for sample selection exists, continue..")
} else {
  cat("Input directory ./Input_Matrix/genes-all_samples-all/ does NOT! exist")
}

dir_out <- paste0(
  dir_base, 
  "Input_Matrix/", 
  "genes-all_", 
  "samples-", 
  ifelse(
    length(remove_samples) > 0,  # remove samples variable is defined in def_variables.R
    paste("no", paste(remove_samples, collapse = "-"), sep = "-"), 
    "all"
  )
)

# create out directory 
if (!file.exists(dir_out)) {
  dir.create(dir_out)
}

# 2.Import files ----
# Check if coldata and metadata are loaded + load if needed

f_load_if_absent_csv(dir_in = dir_in,
                     name_env_variable = "coldata",
                     match_str = "METADATA")

f_load_if_absent_csv(dir_in = dir_in,
                     name_env_variable = "cts",
                     match_str = "CTS")

# 3. Remove undesired samples from dataset ----
# Function to extract row names based on feature_columns and feature_vars
f_samplenames_remove <- function(coldata, feature_columns) {
  
  # Extract columns and features 
  feature_vars <- ls(pattern = "^remove_feature", envir = .GlobalEnv)
  features_rm <- mget(feature_vars, envir = .GlobalEnv)
  names(features_rm) <-feature_columns
  
  matching_rows <- character()  # Initialize an empty vector to collect all matching row names
  
  # Loop through each column in features_rm
  for (feature_col in names(features_rm)) {
    feature_values <- features_rm[[feature_col]]  # Get the corresponding values
    
    # Check if the feature column exists in coldata
    if (feature_col %in% colnames(coldata)) {
      # Find matching row names for each value
      for (value in feature_values) {
        matching_rows <- c(matching_rows, 
                           rownames(coldata[coldata[[feature_col]] == value, ]))
      }
    }
  }
  
  # Remove duplicates and return unique matching rows
  return(unique(matching_rows))
}

list_remove_samples <- f_samplenames_remove(
  coldata = coldata,
  feature_columns = feature_columns
  )

cts_keep <- cts[, !(colnames(cts) %in% list_remove_samples)]
cts_keep <- cts_keep[, order(colnames(cts_keep))] # ensure alphabetic order

coldata_keep <- coldata[!(rownames(coldata) %in% list_remove_samples),]
coldata_keep <- coldata_keep[order(rownames(coldata_keep)),] # ensure alphabetic order


# 3.3 check if the samples in cts and coldata are identical
# Check if the lengths of the vectors are equal
if (
  length(colnames(cts_keep)) != length(rownames(coldata_keep))
) {
  cat("Number of samples in dataframes is unequal.\n")
  cat("Sample removal !NOT! succesful.")
} else {
  # Check if the vectors are identical
  id_cts_coldata <- identical(colnames(cts_keep), rownames(coldata_keep))
  
  if (id_cts_coldata) {
    cat("Sample removal succesful.\n")
    cat("removed", 
        nrow(coldata) - nrow(coldata_keep), 
        "samples for further processing", " \n"
    )
    
  } 
}

# 4. safe files
write.csv(cts_keep, paste0(dir_out, "/", "P1A1_CTS_gen", ".csv"))
write.csv(coldata_keep, paste0(dir_out, "/", "P1A1_METADATA", ".csv"))

cat("Removed from", feature_columns, "\n",
    "Samples matching:", remove_samples, "\n",
    "finished saving files in:", "\n",
    dir_out, "\n")

cts <- cts_keep
coldata <- coldata_keep

rm(list=c("cts_keep", "coldata_keep", "list_remove_samples"))