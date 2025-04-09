f_gen_remove_samples <- function(feature_columns) {
  feature_vars <- ls(pattern = "^remove_feature", envir = .GlobalEnv)
  features_rm <- mget(feature_vars, envir = .GlobalEnv)
  
  # Remove NULL entries from both lists
  valid_idx <- !vapply(features_rm, is.null, logical(1))
  features_rm <- features_rm[valid_idx]
  feature_initials <- substr(feature_columns, 1, 1)[valid_idx]
  
  # If no valid features remain, return empty string
  if (length(features_rm) == 0) {
    return("")
  }
  
  collapsed_tags <- mapply(function(feature_initial, feature_rm) {
    collapsed_features <- paste(feature_rm[lengths(feature_rm) > 1], collapse = "_")
    single_elements <- feature_rm[lengths(feature_rm) == 1]
    combined_features <- c(single_elements, collapsed_features)
    paste0(feature_initial, paste(combined_features, collapse = ""))
  }, feature_initials, features_rm, SIMPLIFY = TRUE)
  
  remove_samples <- paste(collapsed_tags, collapse = "_")
  return(remove_samples)
}


# lists files with a matchable name string and file type in an input directory
f_list_files <- function(dir, match_str, file_type){
  files <- list.files(dir, full.names = F)
  pattern <- paste0(match_str, ".*\\.", file_type, "$")
  file_oi <- files[grep(pattern, files)]
  return(file_oi)
}

# function to read in .csv file in a defined directory, here count and metadata
f_read_csv <- function(dir, csv_file){
  file <- read.table(paste0(dir, csv_file),
                     header = TRUE, 
                     sep = ",", 
                     stringsAsFactors = FALSE, 
                     row.names = 1)
  return(file)
}

# Dynamically generate name based on def_variables 
f_make_name <- function(select_genes, remove_samples) {
  name <- paste0("genes-",
                 if (is.null(select_genes) || length(select_genes) == 0) {
                   "all_"
                 } else {
                   paste0(select_genes, "_")
                 }, 
                 "samples-", 
                 if (remove_samples != "") {  # Check if remove_samples is not an empty string
                   paste("no", remove_samples, sep = "-")
                 } else {  # If remove_samples is an empty string
                   "all"
                 }
  )
  return(name)
}

# f_make_name <- function(select_genes, remove_samples) {
#   name <- paste0("genes-",
#                  if (is.null(select_genes) || length(select_genes) == 0) {
#                    "all_"
#                  } else {
#                    paste0(select_genes, "_")
#                  }, 
#                  "samples-", 
#                  if (length(remove_samples) > 0) {
#                    paste("no", paste(remove_samples, collapse = "-"), sep = "-")
#                  } else {
#                    "all"
#                  }
#   )
#   return(name)
# }


# Check if a variable is present in environment, if not, load it, def for CSV!
f_load_if_absent_csv <- function(dir_in, name_env_variable, match_str) {
  
  # Check if the dataset (name_env_variable) is missing and load it
  if (!exists(name_env_variable, envir = .GlobalEnv)) {
    file_to_load <- f_list_files(dir = dir_in, match_str = match_str, file_type = "csv")
    if (length(file_to_load) > 0) {  # Check if file is found
      data <- f_read_csv(dir = dir_in, csv_file = file_to_load)
      assign(name_env_variable, data, envir = .GlobalEnv)  # Dynamically assign to global environment
      cat(paste0(name_env_variable, " loaded successfully. continue...\n"))
    } else {
      cat(paste0("No ", match_str, " file found.\n"))
    }
  } else {
    cat(paste0(name_env_variable, " already loaded. continue...\n"))
}
}

# define in METADATA/coldata column as factor of interest
f_define_factor <- function(coldata, factor) {
  # Check which factors exist in coldata
  existing_factors <- factor[factor %in% colnames(coldata)]
  
  # Use lapply to convert existing columns to factors
  coldata[existing_factors] <- lapply(coldata[existing_factors], as.factor)
  
  # Return the modified coldata
  return(coldata)
}

# function to calculate counts per million (CPM)
f_CPM <- function(raw_cts){
  samples <- colnames(raw_cts)
  cpm <- matrix(0, nrow = nrow(raw_cts), ncol = ncol(raw_cts),
                dimnames = list(rownames(raw_cts), samples)) 
  
  for (sample in samples) {
    reads_sample <- raw_cts[, sample] # updated from [[]] because its an mtx
  total_sample <- sum(reads_sample)
    cpm_sample <- reads_sample*1/total_sample*10^6  # Formula to calculate CPM
    
    # Store CPM values in the corresponding column of the cpm matrix
    cpm[, sample] <- cpm_sample
  }
  
  return(cpm)
}

# Function to filter input df/matrix by a name genes of interest (goi)
f_filter_df <- function(df, goi_names) {
  df_filter <- df[row.names(df) %in% goi_names, ]
  return(df_filter)
}

# Define the function to save dataframe filtered for goi
f_safe_df <- function(df, name_df, dir, gois, part) {
  file_path <- file.path(dir, 
                         paste0(paste(part, name_df, gois, sep = "_"), 
                                ".csv")
  )
  write.csv(df, file_path, row.names = TRUE)
  return(paste("Saved at", dir))
}

# Load a DEG file and assign name "DEG"
f_assign_load <- function(dir) {
    loaded <- load(dir)
  if (is.character(loaded)) {
    DEG <- get(loaded)
    cat("Data loaded and assigned to variable: DEG")
    return(DEG)
  } else {
    stop("Failed to load data from the specified file.")
  }
}

# Safe ggplots with directory
f_save_plot <- function(dir, name, filetype, plot, width, hight){
  ggsave(paste0(dir,
                name,
                ".",
                filetype
  ),
  plot,
  width = width, # 10
  height = hight # 8
  )
  return <- paste0(paste("saved", name, "at:", sep = " "), dir)
  
  return(return)
} 


# Rename samples to timepoints
f_rename <- function(name) {
  name <- gsub("^P1", "t1", name)
  name <- gsub("^P2", "t2", name)
  name <- gsub("^P3", "t3", name)
  name <- gsub("^P4", "t4", name)
  name <- gsub("^P5", "t5", name)
  name <- gsub("^P6", "t6", name)
  name <- gsub("^P7", "t7", name)
  name <- gsub("^S", "t0", name)
  return(name)
}

# check length object
f_check_length <- function(object, threshold) {
  if (length(object) == 0) {
    message("No matching file!")
  } else if (length(object) == threshold) {
    print("ONE matching file, ready to go....")
  } else if (length(object) > threshold) {
    print("More than one matching file, specify selection!")
  } else {
    print("Some other condition")
  }
}

# read .txt file function
f_read_txt <- function(dir, txt_file){
  file <- read.table(file.path(dir, txt_file),
                     header = TRUE, 
                     sep = "\t",  # Assuming the file is tab-separated. Change this if the separator is different.
                     stringsAsFactors = FALSE, 
                     row.names = 1)
  return(file)
}
