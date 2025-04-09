# Goal:   
# Input:  
# Output: 

# 1.1 clean environment and make sure all graphic devices are closed ----
rm(list = ls())

while (dev.cur() != 1) {
  dev.off()}


# 1.2 Load packages
# library(WGCNA)
# library(tidyverse)
# library(reshape)
library(ggplot2)
library(gridExtra)

# 1.2 Naming
# 1.1 load def_variables.R and def_functions.R
source("/media/wolf2/TKI2019/melanie/Code-transcriptomics-analysis/Pipeline/def_variables.R")
source("/media/wolf2/TKI2019/melanie/Code-transcriptomics-analysis/Pipeline/def_functions.R")

# 1.2 Set wd
setwd(dir_base)



# Define the ranges for count and percent

count <- 1
percent <- 0.5
file <- "_VST_dpi1_count"

# count_range <- c(1, 5, 10)
# percent_range <- c(0.3, 0.4, 0.5, 0.6, 0.7)

# for (count in count_range) { # activate for loop to run over range of counts
#   for (percent in percent_range) {

# 1.3 Directories
# 1.3.1 describe files based on filter criteria
filter <- paste0(file, "_", count, "_", "perc_", percent
                 )

set_name <- f_make_name(select_genes = select_genes,
                        remove_samples = remove_samples
                        )

# 1.3.2 In directory with count matrix
dir_in_counts <- paste0(dir_base, 
                        "Input_filter/",
                        set_name,
                        "/"
                        )


# 1.3.3 In directory for WGCNA
dir_WGCEA <- paste0(dir_base, 
                    "WGCEA/"
                    )

dir_in <- paste0(dir_WGCEA,
                  f_make_name(
                    select_genes = select_genes,
                    remove_samples = remove_samples
                    ),
                  "/"
)

dir_out <- dir_in

# 1.3.4 In genes of interest
dir_goi_base <- paste0(dir_base, 
                  "DGE-goi", 
                  "/", 
                  set_name, 
                  "_", 
                  gois
)

# 2. Load files ----
# 2.1 normalised counts
file_den <- f_list_files(dir = dir_in_counts,
                         match_str = filter,
                         file_type = "csv"
                         )

den <- f_read_csv(dir = dir_in_counts,
                  csv_file = file_den
                  )

# 2.2 load modules WGCEA
name_file <- paste(file, count, "perc", percent, "gene_modules", sep = "_"
                   )

module_df_name <- f_list_files(dir = dir_in,
                               match_str = name_file,
                               file_type = "txt"
                               )

f_check_length(object = module_df_name,
               threshold = 1)

module_df <- f_read_txt(dir = dir_in,
                        txt_file = module_df_name
                        )

# 2.3 load info on goi
goi_info <- read_excel(goi_path)
names(goi_info)[1] <- "goi_ID"
goi_info$goi_ID <- sub("-T1$", "", goi_info$goi_ID)

cat("all files loaded sucesfully...")

# 3. Filter module_df for goi
module_goi <- module_df[module_df$gene_id %in% goi_info$goi_ID,]   

# 3.2 safe the module information
write.table(module_goi,
            file = paste0(dir_out, filter, "_", "effector_modules.txt"),
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

# 4. Short summary
paste("cluster info:")
paste(unique(module_df$colors))
paste("nr effectors in cluster:", nrow(module_goi))

# 5. Visualise results
f_module_expression <- function(modules, module_df, expression_mtx){
  
  submod <- module_df %>%
    subset(colors %in% modules)
  
  row.names(module_df) = module_df$gene_id
  
  # subset expression matrix for modules of interest
  subexpr <- expression_mtx[submod$gene_id,]
  
  return(subexpr)
}

f_plot_module <- function(subexpr, module_df, y_axis_label) {
  
  # Set row names before subsetting
  row.names(module_df) <- module_df$gene_id
  
  # Convert the expression data to a long format for ggplot
  submod_df <- data.frame(subexpr) %>%
    mutate(
      gene_id = row.names(.)
    ) %>%
    pivot_longer(-gene_id) %>%
    mutate(
      module = module_df[gene_id,]$colors
    )
  
  # Plot expression data
  p <- submod_df %>% ggplot(aes(x = name, y = value, group = gene_id)) +
    geom_line(aes(color = module), alpha = 0.2) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) +
    facet_grid(rows = vars(module)) +
    labs(x = "treatment", y = y_axis_label)
  
  return(p)
}

col_oi <- unique(module_df$colors)

for (colour in col_oi){
  
  subexpr <- f_module_expression(modules = colour,
                                 module_df = module_goi,
                                 expression_mtx = den)
  
  p_expr_module <- f_plot_module(subexpr = subexpr,
                                 module_df = module_goi,
                                 y_axis_label = "vst expression")
  
  f_save_plot(dir = dir_out,
              name = paste0("expr_module", colour, filter, "_", gois),
              filetype = "png",
              plot = p_expr_module,
              width = 8,
              hight = 4)
  
  cat("module expression plot saved at:", dir_out)
  
}
