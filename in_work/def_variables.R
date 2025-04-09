# define working directory 
dir_base <- "/home/melanie/working_directory/maize_RNAseq/results/"

# If applicable define which samples should be REMOVED before constructing dds
# ensure the order of feature type and remove samplse matches in order

feature_columns <-c("Treatment", "Time") # "Timepoint" Tissue "Treatment"
remove_feature1 <-c()
remove_feature2 <-c() # "T0", "T1", "C"

# add additional feature rows if needed e.g. remove_featureN, 
# featureN must correspond to feature_columns element N

# Define character (sub)strings to identify which genes derive from which organism
select_genes <- c()

# Define factor(s) of interest for constructing dds
factor_experiment <- c("Treatment", "Time")

# Expression to filter Countdata before DGE
filter_criterion <- "cts[rowSums(cts >=5) >= (20/100)*ncol(cts),]"

# Reference group for DGE/buidling DDS (also relevant for WGCNA = no design)
design_formula <- "~1" # ~1 = No design, ~Treatment
design_shortname <- "NA" # NA = No design, 
ref_gr <- "NA" # NA = no ref gr/design, C

# DGE with DESEeq2: define cutoff lfc, pval for filtering DEG
log2FC <- 1.5 
pval <- 0.05 # padj

# WGCNA parameters
filter_min_count <- 10
filter_min_percent_samples  <- 0.8
BLOCKSIZE <- 5000
POWER <- 12 # 12 = default signed network, 6 = unsigned nework




