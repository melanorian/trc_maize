# This script allows the selection of specific samples from the data set
# The output is a .csv file including the raw counts in the samples of interest

# 0. base set up install/load libraries, setwd, read files ----
# 0.1 load libraries
library(dplyr)
library(ggplot2)
library(PCAtools)
#library("multcompView")
library("purrr")
library(tidyr)
library(ggplot2)
library(DESeq2)

variables_colourd <- c("Treatment", "Time")

mycolors_time <- c(
"#999999","black", "#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2",
"#D55E00", "#CC79A7"
                 )

# Define in and out directories
dir_in <- paste0(
dir_base, "Input_Matrix/",
f_make_name(
  select_genes = select_genes,
  remove_samples = remove_samples
  ),
"/"
)

dir_eda <- paste0(
dir_base, 
"EDA_basic/"
)

name_out <- f_make_name(
select_genes = select_genes,
remove_samples = remove_samples
)

dir_out <- paste0(
dir_eda,
name_out,
"/"
)

if (!file.exists(dir_eda)) {
dir.create(dir_eda)
}

if (!file.exists(dir_out)) {
dir.create(dir_out)
}

f_load_if_absent_csv(dir_in = dir_in,
                   name_env_variable = "coldata",
                   match_str = "METADATA")

f_load_if_absent_csv(dir_in = dir_in,
                   name_env_variable = "cts",
                   match_str = "CTS")


coldata$sample_group <- paste0(coldata$Treatment, "_", coldata$Time)
coldata$Time <- as.factor(coldata$Time)
coldata$Treatment <- as.factor(coldata$Treatment)

dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData = coldata,
  design = as.formula(eval(parse(text = "~ Treatment + Time"))) # design -1 for all variation, with design for specific biologically relevant variation "~ Treatment + Time" or   ~ Treatment + Time + Treatment:Time
)

cts <- assay(vst(dds, blind = FALSE))

# Calculate PCA
pca <- pca(cts, metadata = coldata)

# Make Screeplot
screeplot <- screeplot(pca)

name_plot <- "DesignTT_vst"

plot_name <- paste0(dir_out, "screeplot_", paste0(name_plot), ".svg")

ggsave(plot_name,
       screeplot,
       device = "svg",
       width = 14,
       height = 5)

# Make PCA plot
b <- biplot(pca,
            x = "PC2",
            y = "PC3",
            colby = "Time", 
            colkey = mycolors_time,
            lab = rownames(coldata),
            shape = "Treatment",
            hline = 0, vline = 0,
            legendPosition = 'right',
            showLoadings = FALSE,
            pointSize = 10)


save_plot <- paste0(dir_out, "biplot_", paste0(name_plot), ".svg")

ggsave(save_plot, 
     b,
     device = "svg",
     width = 15,
     height = 15)

return(paste("Saved biplot at",save_plot, "\n")) 


# Safe the created plots 

f_safe_scree(scree_plot = screeplot, 
           name_plot = "CTS_"
           )

# biplot
map2(
.x = pca_plots,
.y = variables_colourd, 
.f = f_safe_biplot,
name_plot = "CTS_"
)

cat("PCA plots saves at:", dir_out)

################################################################################
################################################################################

# 2. import the dds file
file_in <- paste0(dir_in, paste0("DDS_",
                                 f_make_name(select_genes = select_genes,
                                             remove_samples = remove_samples)),
                  ".RData"
)


load(file_in)

# 3. Extracts and order input data ----
# 3.1 Extract raw counts & coldata from dds
cts <- counts(dds)
coldata <- colData(dds)

# 3.2 Use order function to sort based on sample_order in def_variables.R
sorting_index <- order(match(substr(rownames(coldata), 1, 2), sample_order))

# 3.3 Sort based on custom sorting index
coldata_sorted <- coldata[sorting_index, ]
cts_sorted <- cts[, sorting_index]

# 4. PCA plot
# 4.1 calculate PC using raw counts ----
vst <- vst(cts_sorted)

cat("calculated vst...")

count_tables <- list(cts = cts_sorted, vst = vst)
names_count_table <- names(count_tables)

f_pca <- function(count_table, coldata){
  pca <- pca(count_table, metadata = coldata)
  return(pca)
}

pca <- lapply(count_tables, f_pca, coldata = coldata_sorted)

cat("made PCAs...")

# 5. Screeplot
scree_plot <- lapply(pca, screeplot)

f_safe_scree <- function(scree_plot, name_counts){
  ggsave(paste0(dir_out, "screeplot_", paste0(name_counts), ".svg"), 
         scree_plot,
         device = "svg",
         width = 17,
         height = 5)
  return("saved")
}

# Map over both scree_plot and names_count_table simultaneously
scree_plot_saved <- map2(scree_plot, names_count_table, f_safe_scree)
cat("saved screeplots...")

# 6. Biplot
# Create a custom color palette
mycolors_time <- c("#999999", "#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7")
pca_lab <- paste0(coldata_sorted$sample_id) # NULL or paste0(coldata_sorted$sample_id)

f_my_biplot <-function(pca, name_counts){
  
  b <- biplot(pca,
              x = "PC1",
              y = "PC2",
              colby = "time", 
              colkey = mycolors_time,
              lab = pca_lab,
              #shape = 'treatment',
              hline = 0, vline = 0,
              legendPosition = 'right',
              showLoadings = FALSE,
              pointSize = 10)
  
  #b$plot <- b$plot + geom_point(size = 10)  # Adjust the size as needed
  
  ggsave(paste0(dir_out, "biplot_", paste0(name_counts, "_", filter), ".svg"), 
         b,
         device = "svg",
         width = 15,
         height = 15)
  
  return("saved biplot")
  
}

# run the biplot function
map2(pca, names_count_table, f_my_biplot)

cat("saved plots in:", dir_out)
