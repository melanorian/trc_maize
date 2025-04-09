# Goal: This script produces a list of differentially expressed genes in the target
# based on defined reference levels 
# Input: Requires a count matrix + metadata file (), dds with selected samples (P1B)
# Output: dds with defined reference group, results of DESeq() for all genes, 
#         results of DESeq() for all genes given filtering cirteria for log2FC and padj,
#         results of DESeq() log2FC/padj filter + filtered for organism of interest


# 1.1 load libraries
library(dplyr)
library(DESeq2)

# 1.3 in direcory
dir_in <- paste0(dir_base, "DDS/")

# 1.4 out directory
dir_dge <- paste0(dir_base, "DGE/")

name_out <- paste0(ref_gr, "_", 
                   f_make_name(select_genes = select_genes,
                               remove_samples = remove_samples) # var defined in main.R
)

dir_out <- paste0(dir_dge, name_out, "/")

if (!file.exists(dir_dge)) {
  dir.create(dir_dge)
}

if (!file.exists(dir_out)) {
  dir.create(dir_out)
}

# 2. import the dds file
file_in <- paste0(dir_in, 
                  paste0(
                    dir_out, "DDS_", 
                    set_name, "_",
                    "Des", design_shortname,
                    ".RData")
)

if (!exists("dds", envir = .GlobalEnv)) {
  load(file_in)
  cat("dds loaded. continue compute DEGs...")
  }else{
    cat("dds loaded. continue compute DEGs...")
}

# 3. Identify differentially expressed genes ----
# Info DGE analysis: https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/05b_wald_test_results.html
# 3.1 Run DESe1 function
dds <- DESeq(dds)

list_contr <- as.vector(resultsNames(dds)) 
cat("Contrasts in DGE:", "\n",
    list_contr) # Check list of contrasts computed

# 4. Extract results per contrast 
res_list <- lapply(list_contr, function(contrast) {
  # Extract results for the current contrast
  result <- results(dds, name = contrast, alpha = pval, pAdjustMethod = "BH")
  return(result)
})

names(res_list) <- list_contr # resulting list comprises all genes
list_contr <- list_contr[-1] # skip intercept

# # 6. Shrink logFC towards 0 to account for low count/high dispersion ----
# res_list <- lapply(list_contr, function(contrast) {
#   # Extract result corresponding to the contrast
#   print(paste("CONTRAST", contrast))
#   result <- res_list[[contrast]]
#   print(paste("RESULT-reslist res happened"))
#   # Apply lfcShrink
#   shrink_result <- lfcShrink(dds, coef = contrast, res = result, type = "apeglm") # SELECT TYPE, only apeglm yields resonable results
#   print(paste("SHRUNK - happened"))
#   return(shrink_result)
# })

names(res_list) <- list_contr

# 3.3 Filter for genes with log2FC, padj within range of interest $ sort in DEGs, UP, Down regulated 
# 3.3.1 write function
f_DEG_list <- function(results_list, contrast) {
  
  # Get the results table for the contrast of interest
  res <- results_list[[contrast]]
  # Extract DEGs based on padj and log2 fold change thresholds
  DEGs <- subset(res, padj < pval & abs(log2FoldChange) >= log2FC)
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

# 3.3.2 Run the function over the list of different inter-time point comparisons
DEG_list <- lapply(list_contr, f_DEG_list, results_list = res_list)
names(DEG_list) <- list_contr

# 4. Filter for organism of interest
# 4.1 Vector storing the different (sub)sets of DEGs in DEG_list
type_DEG <- names(DEG_list[[1]]) # access the names of the element stored in the first element of the DEG list

# 4.2 Define a function to generate a list of DEGs filtered for the organism of interest
# the function takes as dynamic input a list of contrasts and types of DEGs (DEG/up/down)

filtered_DEGs <- lapply(names(DEG_list), function(contrast, DEG_list, type_DEG, str_match) {
  filtered_contrast <- list()
  for (type in type_DEG) {
    DEG_names <- rownames(DEG_list[[contrast]][[type]])
    
    # Filter for gene names that start with the specified substring
    matched_genes <- grep(paste0("^", str_match), DEG_names, value = TRUE)
    
    filtered_DEGs <- DEG_list[[contrast]][[type]][matched_genes, ]
    filtered_contrast[[type]] <- filtered_DEGs
  }
  return(filtered_contrast)
}, DEG_list = DEG_list, type_DEG = type_DEG, str_match = genes) #str_match = genes defined in the begining of the script

# Rename filered_DEG list by comparisons
filtered_DEGs <- setNames(filtered_DEGs, names(DEG_list))
print(names(filtered_DEGs))

# 5. Overview DGE 
# 5.1 Create a matrix to store the counts
count_matrix <- matrix(NA, nrow = length(type_DEG), ncol = length(list_contr), dimnames = list(type_DEG, list_contr))

# 5.2 Fill the matrix with counts
for (contrast in list_contr) {
  for (type in type_DEG) {
    count_matrix[type, contrast] <- nrow(filtered_DEGs[[contrast]][[type]])
  }
}

# 5.3 Convert matrix to data frame
DEG_summary <- as.data.frame(count_matrix)
print(DEG_summary)

# 6. Save output
# 6.1 safe updated dds = dds_out
save_res_list <- paste(
  "P2A", "unfiltered", 
  paste0("FC", "default"), 
  paste0("PADJ", pval),
  paste0("Des", design_shortname),
  sep = "_")

save(res_list, 
     file = paste0(
       dir_out, "/", 
       save_res_list, 
       "Des", 
       design_shortname,
       ".RData")
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

# 6.3 Results of DEG_list filtered for genes only from organism of interest e.g. Pe - filtered_DEGs
save_filtered_DEGs <- paste("P2A", "filterLFC", genes, paste0("FC", log2FC), paste0("PADJ", pval), sep = "_")
save(filtered_DEGs, file = paste0(dir_out, "/", save_filtered_DEGs, ".RData"))

# 6.4 
save_DEGoverview <- paste("P2A", "filter-overview", genes, paste0("FC", log2FC), paste0("PADJ", pval), sep = "_")
write.csv(DEG_summary, file = paste0(dir_out, "/", save_DEGoverview, ".csv"), row.names = TRUE)

# 6.5 Overwrite DDS from Part1B1
save(dds, 
     file = paste0(
       dir_out, "DDS_", 
       set_name, "_",
       "Des", design_shortname,
       ".RData")
)
     

cat("saved DESEq results in: ", dir_out)

rm("")
        