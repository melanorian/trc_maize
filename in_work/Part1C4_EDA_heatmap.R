# 5. Hierarchical clustering ----

library(pheatmap)
library(dplyr)
library(ggplot2)
library("svglite")

# Change Variables -------------------------------------------------------------

# Name matrix to plot if != cts
cts <- cts
name_mtx_in <- "CTS" # name for saving output

# Define Variables for annotation
variables_colourd_discrete <- c("Timepoint", "Treatment", "Tissue")
variables_colourd_continous <- c() #"reads"

# NR of clusters for heatmap
nr_clusters <- 1

# Set colours for annotation of heatmap
col_factor1 <- c(
"#999999", "black", "#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", 
"#D55E00", "#CC79A7", "red", "red", "red", "red"
)

# for factor treatment
col_factor2 <- c(
"#332288", "#CC6677", "red", "red", "red", "red"
)

col_factor3 <- c(
"#DDCC77", "#117733", "red", "red", "red", "red"
)

# Define color palettes for continuous variables
col_contin1 <- colorRampPalette(c(
"#FFCF65", "#E69F00", "#8A5F00", "#362500")
)(1000)

col_contin2 <- colorRampPalette(c(
"#A7FFE8", "#009E73", "#00644A", "#004C38")
)(1000)

col_contin3 <- colorRampPalette(c(
"#E8C2D7", "#CC79A7", "#883462", "#3A162A")
)(1000)

#-------------------------------------------------------------------------------

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

# Check if dir_out already exist, if not make them
if (!file.exists(dir_eda)) {
dir.create(dir_eda)
}

if (!file.exists(dir_out)) {
dir.create(dir_out)
}

# Check if data already present or load
f_load_if_absent_csv(dir_in = dir_in,
                   name_env_variable = "coldata",
                   match_str = "METADATA")

f_load_if_absent_csv(dir_in = dir_in,
                   name_env_variable = "cts",
                 match_str = "CTS")


# Calculate correlation matrix
if (!exists("cor_mtc")) {
# Calculate the Spearman correlation if 'cor' does not exist
cor_mtc <- cor(cts, method = "spearman")
} else {
print("Correlation matrix exists.")
}

# Make df for annotation
coldata_hm <- coldata[, paste(c(
variables_colourd_discrete, 
variables_colourd_continous)
)]

coldata_hm <- coldata_hm %>% 
mutate(across(where(is.character), as.factor))

# To annotate heatmap we have to create a list which assigns colours to the
# factors we want to annotate on the heatmap.
# Check if enough colours available for variable, correct naming of lists 
# with variable of interest. Pretty complainy function..

# Create empty list
mycolors <- list()

# Assign colors for discrete variables to empty list
for (i in seq_along(variables_colourd_discrete)) {
col_factor_name <- paste0("col_factor", i)
factor_name <- colnames(coldata_hm)[i]  # Get the column name for the factor

factor_levels <- levels(coldata_hm[[factor_name]])  # Get the levels of the factor
mycolors[[factor_name]] <- setNames(
  get(col_factor_name)
  [1:length(factor_levels)],
  factor_levels
)}

# Remove Non-assigned colours
mycolors <- lapply(mycolors, function(x) {
x[!is.na(names(x))] 
})

# Add assigned colors for continuous variables to list
for (i in seq_along(variables_colourd_continous)) {
col_factor_name <- paste0("col_contin", i)

start_position <- which(
  !colnames(coldata_hm) %in% 
    variables_colourd_discrete)[1]

continuous_factor_name <- colnames(coldata_hm)[start_position]
mycolors[[continuous_factor_name]] <- get(col_factor_name)
}


# 5.3.2 change annotation colour to a distinguishable scale 
title_heatmap <-  paste(
"Heatmap_cor", 
paste0("clust", nr_clusters), 
sep = "_"
)

# 5.3.3 draw heatmap
h1 <- pheatmap(
cor_mtc, 
annotation = coldata_hm,
annotation_colors = mycolors, 
cluster_cols = TRUE,
cutree_cols = nr_clusters,
border_color = "NA",
main = title_heatmap,
show_rownames = FALSE,
show_colnames = FALSE
)

title_heatmap <-  paste(
"Heatmap_cor", 
paste0("clust", nr_clusters), 
sep = "_"
)

ggsave(paste0(
dir_out, 
"1PC4_", title_heatmap, 
".svg"
),
h1,
device = "svg",
width = 15,
height = 13
)


cat("Hetmap EDA saved at:", dir_out)
