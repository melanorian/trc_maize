#!/usr/bin/env Rscript

# Goal: construct a gene-level countmatrix for further analysis and METADATA file with 
# information on reads/sample
# Input: working directory defined in def-variables containing count_table dir with
# CTS_gene.csv and METADATA.txt
# Output: "Input_Matrix" directory with subdirectory "genes-all_samples-all". Files 
# CTS-trc.csv and METADATA.csv will be generated. Both sorted alphabetically by sample names

dir_in <- paste0(dir_base, "count_table")
dir_out <- paste0(dir_base, "Input_Matrix") 

dir_out_file <- paste0(dir_out, 
                     "/", 
                     "genes-all_samples-all" 
                     ) 

# 1.4 create directories if they don't exist yet
if (!file.exists(dir_out)) {
dir.create(dir_out)
}

if (!file.exists(dir_out_file)) {
dir.create(dir_out_file)
}

# Check if the count directory exists
if (dir.exists(dir_in)) {
cat("Count Directory exists: ", dir_in, "\n")
} else {
cat("Cound Directory does !NOT! exist: ", dir_in, "\n")
}

# 2. load CTS file
file_cts <- paste(dir_in, "CTS_gene.csv", sep = "/")

# Check if count file exists
if (file.exists(file_cts)) {
cat("Count File exists: ", file_cts, "\n")
} else {
cat("Count File does !NOT! exist: ", file_cts, "\n")
}

cts <- read.table(file_cts, 
                header = TRUE, 
                sep = ",", 
                stringsAsFactors = FALSE, 
                row.names = 1
                )

cat("loaded count file...", "\n")

# load METADATA file
file_col <- paste0(dir_in, "/", "METADATA.csv")

# Check if METADATA file exists
if (file.exists(file_col)) {
cat("METADATA exists: ", file_col, "\n")
} else {
cat("METADATA does !NOT! exist: ", file_col, "\n")
}

coldata <- read.table(file_col, 
                header = TRUE, 
                sep = ",", 
                stringsAsFactors = FALSE
)

cat("loaded metadata file...", "\n")

# Load metadata on total nr of reads
file_reads <- paste0(dir_in, "/", "read_nr.csv")

# Check if READ file exists
if (file.exists(file_col)) {
cat("file with nr of reads per sample exists: ", file_col, "\n")
} else {
cat("file with nr of reads per sample does !NOT! exist: ", file_col, "\n")
}

reads <- read.table(file_reads,
                  header = TRUE,
                  sep = ",",
                  stringsAsFactors = FALSE
)

cat("file with nr of reads per sample...", "\n")

# 3. Summarize to gene level with multiple transcripts
# 3.1 Use common GeneID to summarize to gene level
gennames <-  gsub("_.*", "", rownames(cts)) # Extract geneID

cat("transcript names head:",
head(rownames(cts)), 
sep = " ",
"\n"
)

cat("gene names head:",
head(gennames), 
sep = " ",
"\n"
)

# Sum different transcripts to gene level by geneID
cts_genlevel <- aggregate(. ~ gennames, data = cts, FUN = sum)
row.names(cts_genlevel) <- cts_genlevel$gennames
cts_genlevel <- cts_genlevel[,-1]
cts_genlevel <- cts_genlevel[, order(colnames(cts_genlevel)), drop = FALSE]
colnames(cts_genlevel) <- gsub("^X", "", colnames(cts))

# 4. METADATA
rownames(coldata) <- coldata$Sample
coldata <- coldata[order(row.names(coldata)), , drop = FALSE]
samplename <- rownames(coldata)

# 4.2 add column read nr extracted from fastq files
rownames(reads) <- reads$sample
reads <- reads[order(row.names(reads)), , drop = FALSE]

# 3.2 Print outputs
cat("column names and order CTS", "\n",
  colnames(cts_genlevel), "\n"
  ) 

cat("nr multi transcripts:", 
  nrow(cts) - nrow(cts_genlevel), 
  "\n"
  ) # check nr rows pre-summarizing

# Rename df by Group for easier understanding of output
coldata$ID <- paste0(coldata$Treatment, "_", sprintf("%.1f", coldata$Time_h), "_", coldata$Biol_repl) # Create the ID column
coldata$ID <- gsub("\\.", "", coldata$ID)
rownames(coldata) <- coldata$ID

# use the same names for cts df
colnames(cts_genlevel) <- rownames(coldata)

# names for reads
rownames(reads) <- rownames(coldata)

# merge mapping info from reads file to metadata
coldata <- full_join(coldata, reads)
rownames(coldata) <- coldata$ID

# Check if the lengths of the vectors are equal
if (
length(colnames(cts_genlevel)) != length(rownames(reads))| 
length(colnames(cts_genlevel)) != length(rownames(coldata))
) {
cat("Number of samples in dataframes is unequal!.\n")
} else {
# Check if the vectors are identical
id_cts_read_nr <- identical(colnames(cts_genlevel), rownames(reads))
id_cts_coldata <- identical(colnames(cts_genlevel), rownames(coldata))
id_read_nr_coldata <- identical(rownames(reads), rownames(coldata))

if (id_cts_read_nr && id_read_nr_coldata) {
  cat("All vectors are identical.\n")
} else {
  if (!id_cts_read_nr) {
    cat("Number samples in: count data != read counts.\n")
  }
  if (!id_cts_coldata) {
    cat("Number samples in: count data != metadata\n")
  }
  if (!id_read_nr_coldata) {
    cat("Number samples in: metadata != read counts .\n")
  }
}
}

cat("safe output...", "\n")

# Round counts
cts_genlevel <- round(cts_genlevel)

# 5. Safe Gene level cts output for further work
write.csv(cts_genlevel, paste0(dir_out_file, "/", "P0_", "CTS_gen", ".csv"))
write.csv(coldata, paste0(dir_out_file, "/", "P0_", "METADATA", ".csv"))

cat("as:", paste0(dir_out_file, "/", "P0_", "CTS_gen", ".csv"), "\n")
cat("and:", paste0(dir_out_file, "/", "P0_", "METADATA", ".csv"), "\n")

cts <- cts_genlevel
rm(list = c("reads", "cts_genlevel", "gennames"))
