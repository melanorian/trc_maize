# Goal: This script produces a list of differentially expressed genes in the target
# based on defined reference levels 
# Input: Requires a count matrix + metadata file (), dds with selected samples (P1B)
# Output: dds with defined reference group, results of DESeq() for all genes, 
#         results of DESeq() for all genes given filtering cirteria for log2FC and padj,
#         results of DESeq() log2FC/padj filter + filtered for organism of interest


# 1.1 load libraries
library(dplyr)
library(DESeq2)
library(purrr)

# SET THESE VARIABLES
# List all timepoints for analysis
days <- c("T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8")
design_formula <- as.formula(~"Treatment")

# # This greates based
# group <- c(paste0(
#   "t",
#   as.numeric(gsub(".*?(\\d)$", "\\1", ref_gr)) + 1
# ))

# 1.3 in direcory
dir_in_raw <- paste0(dir_base, "Input_Matrix/",
                     f_make_name(select_genes = select_genes,
                                 remove_samples = remove_samples),
                     "/"
)
# 1.4 Create all Output directories/ check if they already exist

name_out <- paste0(ref_gr, "_", 
                   f_make_name(select_genes = select_genes,
                               remove_samples = remove_samples)
)

dir_dge <- paste0(dir_base, "DGE/")
dir_dds <- paste0(dir_base, "DDS/")
dir_dge_consecutive_comps <- paste0(dir_dge, "consecutive_comps/")
dir_dds_consecutive_comps <- paste0(dir_dds, "consecutive_comps/")
dir_out_deg <- paste0(dir_dge_consecutive_comps, name_out, "/")
dir_out_dds <- paste0(dir_dds_consecutive_comps, name_out, "/")

if (!file.exists(dir_dge)) {
  dir.create(dir_dge)
}

if (!file.exists(dir_dds)) {
  dir.create(dir_dds)
}

if (!file.exists(dir_dge_consecutive_comps)) {
  dir.create(dir_dge_consecutive_comps)
}

if (!file.exists(dir_dds_consecutive_comps)) {
  dir.create(dir_dds_consecutive_comps)
}

if (!file.exists(dir_out_deg)) {
  dir.create(dir_out_deg)
}

if (!file.exists(dir_out_dds)) {
  dir.create(dir_out_dds)
}
# 2.Import files ----
# Check if coldata and metadata are loaded + load if needed

f_load_if_absent_csv(dir_in = dir_in_raw,
                     name_env_variable = "coldata",
                     match_str = "METADATA")

f_load_if_absent_csv(dir_in = dir_in_raw,
                     name_env_variable = "cts",
                     match_str = "CTS")

cat("loaded input data from:", dir_in_raw)

# 3. Prepare data for constructing DESeq2 object ----
# 3.1 Ensure that columns of interest for DGE are factors (for DESeq())
coldata <-  f_define_factor(coldata, factor_experiment[1])

cat("Levels for DGE of consecutive timepoints:",
    levels(coldata[,factor_experiment[1]])
    )

# 3.2 Construct dds dynamically
# 3.2.1 Function for making individual dds for each consecutive variable
f_dds_consecutive <- function(day, coldata, cts){
  # Step1: Filter METADATA, cts by day of interest
  coldata_keep <- coldata[coldata$Timepoint == day, ]
  cts_keep <- cts[, coldata$Timepoint == day]
  # Step2: Creat DDS 
  dds <- DESeqDataSetFromMatrix(
    countData = cts_keep,
    colData = coldata_keep,
    design = ~Treatment
    )
  
  return(dds)
}

# 3.2.2 Apply function to create list of dds 
dds_list <- lapply(
  days, 
  f_dds_consecutive,
  coldata = coldata,
  cts = cts
  )

# 3.2.3 Name elements in dds by matching consecutive variable 
names(dds_list) <- days

# 4.Run DGE analysis - THIS TAKES TIME! Run only if needed!

# 4.1 Function to run DESeq and assign results back to dds
f_deseq_consecutive <- function(dds){
  dds <- DESeq(dds)
  return(dds)
}

# 4.2  Apply DESeq() function to conduct DGE analysis. 
deseq_list <- lapply(dds_list, f_deseq_consecutive)

# 4.3 Rename deseq_list elements by consecutive variable
names(deseq_list) <- days

# 4.4. List contrasts stored in results, ignore intercept = first element
list_contr <- as.vector(resultsNames(deseq_list[[1]]))[-1]

# 5. Extract information on DEGs from the DEseq results ----

# 5.1 Function to extract results from dds (result of DESeq())
f_res_consecutive <- function(deseq_dds, contrast){
  # Step1: Extract results for the current contrast
  # Variables defined in def_variables for consistency
  result <- results(deseq_dds, 
                    name = contrast, 
                    alpha = pval, 
                    pAdjustMethod = "BH"
                    )
  return(result)
}

# 5.2 Apply function to extract results for each DESeq() dds
res_list <- lapply(deseq_list, 
                   f_res_consecutive,
                   contrast = list_contr # (here only length 1, >1 need map2)
                   )

# 6. Sort information from results of DGE for further analysis
# 6.1 Function to extract results larger than defined logFC, with specific pval
#     split in up/downregulated genes
f_DEG_list <- function(result_in, pval, log2FC) {
  # Extract DEGs based on padj and log2 fold change thresholds
  DEGs <- subset(result_in, padj < pval & abs(log2FoldChange) >= log2FC)
  # Split DEGs into upregulated and downregulated genes
  upregulated_genes <- subset(DEGs, log2FoldChange >= log2FC)
  downregulated_genes <- subset(DEGs, log2FoldChange <= -log2FC)
  # Make a list of up/downregulated genes
  DEGs_list <- list(DEGs, upregulated_genes, downregulated_genes)
  # Rename the list elements
  names(DEGs_list) <- c("DEGs", "up", "down")
  # Return the DEGs list
  return(DEGs_list)
}

# 6.2 Run the function over the list of different inter-time point comparisons
DEG_list <- lapply(res_list, 
                   f_DEG_list, 
                   pval = pval,
                   log2FC = log2FC
                   )

# 7. Quick overview over results ----
# 7.1 Create a matrix to store the counts
summary_matrix_DEG <- matrix(
  NA, 
  nrow = length(type_DEG), 
  ncol = length(days), 
  dimnames = list(type_DEG, days)
  )

# 5.2 Fill the matrix with counts
for (contrast in days) {
  for (type in type_DEG) {
    summary_matrix_DEG[type, contrast] <- nrow(DEG_list[[contrast]][[type]])
  }
}

# 7.3 Convert matrix to data frame
DEG_summary <- as.data.frame(summary_matrix_DEG)
print(DEG_summary)


# 8. Save output ---
# 8.1 Safe dds with results of DESeq2
# shared name for all results
# save_res_list <- paste(
#   "P2A", "unfiltered", 
#   paste0("FC", "default"), 
#   paste0("PADJ", pval),
#   paste0("Design", design_shortname),
#   sep = "_")

# Funciton to safe files seperately seperately as RData
f_safe_dds_consecutive <- function(results_in, day, dir_out, name_file){
  
  file_name <- paste0(
    dir_out,
    name_file,
    "_",
    day,
    ".RData")
  
  save(results_in, file = file_name)
  
  return(file_name)
  }

# Safe individual dds
map2(.x = deseq_list,
     .y = days,
     .f = ~f_safe_dds_consecutive(
       results_in = .x,
       day = .y,
       dir_out = dir_out_dds,
       name_file = "DDS"
     )
     )

# Safe unfiltered res list sperately
map2(.x = res_list,
     .y = days,
     .f = ~f_safe_dds_consecutive(
       results_in = .x,
       day = .y,
       dir_out = dir_out_deg,
       name_file = "DGE_RESULTS"
     )
)

# Safe unfiltered res list joined
save(res_list, 
     file =paste0(dir_out_deg,
                  "DGE_RESULTS-joined",
                  ".RData"
     )
)

# Safe List of DEGs as one Rdata file
save(DEG_list, 
     file =paste0(dir_out_deg,
                  paste("DGE_RESULTS-joined-filterd",
                        paste("padj", pval, sep = "-"),
                        paste("logFC", log2FC, sep = "-"),
                        sep = "_"
                        ),
                  ".RData"
                  )
     )

# 6.2 Results for DEG based on DEFINED parameters of log2FC, padj - DEG_list
save_DEG_list <- paste("P2A", 
                       "filterLFC", 
                       paste0("FC", log2FC), 
                       paste0("PADJ", pval), 
                       paste0("Des", design_shortname),
                       sep = "_"
                       )

save(DEG_list, file = paste0(dir_out, 
                             "/", 
                             save_DEG_list, 
                             "_" ,"Des", 
                             design_shortname,".RData")
     )

# 6.4 Save overview 
write.csv(DEG_summary, 
          file = paste0(
            paste0(
              dir_out_deg,
              paste(
                "Overview", 
                paste0("FC", log2FC), 
                paste0("PADJ", pval), 
                sep = "_"
                ),
              ".csv"
              ) 
            ), row.names = TRUE
          )

cat("saved DESEq results over consecutive variable in: ", "\n",
    dir_out_dds, "\n",
    "AND", "\n",
    dir_out_deg, "\n"
)