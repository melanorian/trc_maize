# 1.1 clean env
rm(list = ls())

while (dev.cur() != 1) {
dev.off()}

dir_code <- "/home/melanie/working_directory/maize_RNAseq/trc_maize/in_work/"

# 1.2 Source shared variables from def_variables.R script
source(paste0(dir_code ,"def_variables.R"))

# 1.3 set working directory
setwd(dir_base) # def in def_variables.R
cat("Working directory: ", "\n", getwd(), "\n") # print current wd to terminal

# 1.4 load shared functions
source(paste0(dir_code ,"def_functions.R"))

# generate tag for removing samples if defined
remove_samples <- f_gen_remove_samples(feature_columns)

# 2. Pre-processing of data ----

# 2.1 Check if directory and files already exist
dir_all_check<- paste0(dir_base, "Input_Matrix/genes-all_samples-all")
all_cts_check <- paste0(dir_all_check, "/P0_CTS_gen.csv")
all_metadata_check <- paste0(dir_all_check, "/P0_METADATA.csv")

# Check if the directory and both files exist
if (!dir.exists(dir_all_check) || !file.exists(all_cts_check) || !file.exists(all_metadata_check)) {

cat("Convert transcript level to gene level counts... \n")

# Process transcript level to gene level counts
source(paste0(dir_code ,"Part0_preprocessing.R"))
} else {
cat("Transcript level count tables and metadata file already exist.\n")
}

# 3. Sample/Gene selection ----

# Only run this section if samples and/or genes need to be removed
# It looks only complicated because I didn't wanna hard code too many things so
# that they can easily be changed from def_variables.R. Also, I don't want to run
# computationally intense secteions if the output already exists. That's why all
# the if statements are taking up so many lines ;)

# 3.1 Directory/file where output of sample/gene selection is saved
dir_out_Part1A <- paste0(
  dir_base, "Input_Matrix/", "genes-all_", "samples-", 
  ifelse(
    nzchar(remove_samples),  # remove samples variable is defined in def_variables.R
    paste("no", paste(remove_samples, collapse = "-"), sep = "-"), 
    "all"
  )
)

file_out_Part1A <- paste0(
  dir_out_Part1A, "/",
  f_list_files(
    dir = dir_out_Part1A,
    match_str = "CTS",
    file_type = "csv"
  )
)

# 3.1 Check very lengthy if the output files already exist only run the scripts
# if the output does not exist yet.

if (length(feature_columns) == 0 && length(select_genes) == 0) {
  
  cat("Skipping all sample selection \n")
  
} else if (length(feature_columns) > 0 && length(select_genes) == 0) {
  
  if (!file.exists(file_out_Part1A)) {
    
    cat("Remove", remove_samples, "samples", "\n")
    source(paste0(dir_code, "Part1A1_select_samples.R"))
    
  } else {
    
    cat("File with removed samples already exists at:", "\n",
        file_out_Part1A, "\n",
        "Skipping Part1A1...", "\n"
    ) 
  }
  
} else if (length(select_genes) > 0 && length(feature_columns) == 0) {
  
  if (!file.exists(file_out_Part1A)) {
    
    cat("Select", select_genes, "for further analysis..\n")
    source(paste0(dir_code, "Part1A2_select_samples.R"))
    
  } else {
    
    cat("File", file_out_Part1A, "already exists. Skipping Part1A2.\n")
  }
  
} else if (length(feature_columns) > 0 && length(select_genes) > 0) {
  
  if (!file.exists(file_out_Part1A)) {
    
    cat("Remove samples:", remove_samples, "\n")
    source(paste0(dir_code, "Part1A1_select_samples.R"))
    
  } else {
    
    cat("File", file_out_Part1A, "already exists. Skipping Part1A1.\n")
  }
  
  if (!file.exists(file_out_Part1A)) {
    
    cat("Select genes:", select_genes, "\n")
    source(paste0(dir_code, "Part1A2_select_samples.R"))
    
  } else {
    
    cat("File", file_out_Part1A, "already exists. Skipping Part1A2.\n")
  }
}

# # 4 EDA optional - Too individual per dataset/time intense too run from here
#     check out filtering for low read counts
# #source(paste0(dir_code, "Part1C3_EDA_PCA.R"))
# source(paste0(dir_code, "Part1C4_EDA_heatmap.R"))
# source(paste0(dir_code, "Part1C5_EDA_genes_heatmap.R"))

# 5.Construct DDS, extract normalised counts, filtered count matrix

# check if initial directory already exist matrix with specifi combi of 
# filtering to check if the script was run P1B1_cts_filtered_5_20.csv
# Extract cutoff and percentage of samples before the if condition
###-----

if (!file.exists(
  
  # Construct output filepath
  paste0(
    dir_base, 
    "DDS/",
    "DDS_", 
    f_make_name(select_genes = select_genes, 
                remove_samples = remove_samples), "_",
    "Des", design_shortname,
    ".RData")
  
)) {
  # Create a local block to avoid global environment pollution
  {
    # Run the source script
    cat("Running the script to generate:", "\n")
    source(paste0(dir_code, "Part1B1_DGE_makeDESeq2_dataset.R"))
  }
  
} else {
  cat("...Output file already exists here:", "\n", 
      
      paste0(
        dir_base, 
        "DDS/",
        "DDS_", 
        f_make_name(select_genes = select_genes, 
                    remove_samples = remove_samples), "_",
        "Des", design_shortname,
        ".RData"),
      "\n",
      
      "...Skipping script: Part1B1_DGE_makeDESeq2_dataset.R.... \n")
}


# 6. Run DESeq2 analysis
if (!file.exists(
  
  # Construct output filepath
  paste0(
    dir_base, "DGE/",
    paste0(ref_gr, "_", 
           f_make_name(
             select_genes = select_genes,
             remove_samples = remove_samples),
           "/",
    paste("P2A", "filter-overview", 
          genes, paste0(
            "FC", log2FC), 
          paste0(
            "PADJ", pval),
          sep = "_"
          ),
    ".csv"
    )
  )
)) {
  # Create a local block to avoid global environment pollution
  {
    # Run the source script
    cat("Running DESeq2 to identify DEGs:", "\n")
    source(paste0(dir_code, "Part2A_DGE_compute.R"))
  }
  
} else {
  cat("DESeq was already run for this input.Skipping script.\n")
}



# 7. WGCNA
# 7.1 Calculate optimal soft thrashold for WGCNA 

if(file.exists(paste0(
  paste0(dir_WGCEA,
         f_make_name(
           select_genes = select_genes,
           remove_samples = remove_samples),
         "/",
         filter_min_count,
         filter_min_percent_samples,
         "grid_power",
         ".png"
         )
  )
  )

){
  cat("Plot to choose soft thrashold parameter for WGCNA exists.", "\n",
      "check at:", paste0(dir_WGCEA,
                          f_make_name(
                            select_genes = select_genes,
                            remove_samples = remove_samples),
                          "/",
                          filter_min_count,
                          filter_min_percent_samples,
                          "grid_power",
                          ".png"
      ))
}else{
  cat("calculating optimal soft thrashold parameter for WGCNA", "\n")
  source("/home/melanie/potato_RNAseq/RNAseq_potato_waterlogging/DGE_pipeline/Part5A_WGCEA_explore_modelfit.R")
}


# 7.2 Run WGCNA 

if(file.exists(
  paste0(
    dir_base,
    "WGCEA/",
    f_make_name(select_genes = select_genes,
            remove_samples = remove_samples
            ),
    "/",
    "WGCNA_RES_", 
    "count", filter_min_count, 
    "perc", filter_min_percent_samples, 
    "_power", POWER, 
    "_block", BLOCKSIZE, 
    ".RDS"
  ))){
  cat("WGCNA output already exists at:", "\n",
      paste0(
        dir_base,
        "WGCEA/",
        f_make_name(select_genes = select_genes,
                    remove_samples = remove_samples
        ),
        "/",
        "WGCNA_RES_", 
        "count", filter_min_count, 
        "perc", filter_min_percent_samples, 
        "_power", POWER, 
        "_block", BLOCKSIZE, 
        ".RDS"
      ))
}else{
  cat("Run WGCNA analysis and safe output")
  source("/home/melanie/potato_RNAseq/RNAseq_potato_waterlogging/DGE_pipeline/Part5B_run_WGCEA.R")
}


if(file.exists(
  paste0(
    dir_base,
    "WGCEA/",
    f_make_name(select_genes = select_genes,
                remove_samples = remove_samples
                ),
    "EIGENGENES_", 
    "count", filter_min_count, 
    "perc", filter_min_percent_samples, 
    "_power", POWER, 
    "_block", BLOCKSIZE,  
    "_GRID.png"
)
  
)){
  cat("Eigengen plots of WGCNA alraedy exist at:", "\n",
      paste0(
        dir_base,
        "WGCEA/",
        f_make_name(select_genes = select_genes,
                    remove_samples = remove_samples
        ),
        "EIGENGENES_", 
        "count", filter_min_count, 
        "perc", filter_min_percent_samples, 
        "_power", POWER, 
        "_block", BLOCKSIZE,  
        "_GRID.png"
      ))
}else{
  cat("Plot Eigengenes from WGCNA. Continue sourcing script....", "\n")
  source("/home/melanie/potato_RNAseq/RNAseq_potato_waterlogging/DGE_pipeline/Part5C_WGCEA_eigengene.R")
}


