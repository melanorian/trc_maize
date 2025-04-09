# This script allows the selection of specific samples from the data set
# The output is a .csv file including the raw counts in the samples of interest
# 0. base set up install/load libraries, setwd, read files ----
# 0.1 load libraries
library(tidyr)
library(ggplot2)

# 1.3 in direcory
#dir_in <- paste0(dir_base, "DDS/")

dir_in <- paste0(dir_base, "Input_Matrix/",
                 f_make_name(
                   select_genes = select_genes,
                   remove_samples = remove_samples)
                 ,
                 "/"
                 )

# 1.4 out directory
name_out <- f_make_name(
  select_genes = select_genes,
  remove_samples = remove_samples
  )

dir_eda <- paste0(
  dir_base, 
  "EDA_basic/"
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


# 2. import the dds file
# file_in <- paste0(dir_in, paste0("DDS_",
#                  f_make_name(select_genes = select_genes,
#                              remove_samples = remove_samples)),
#                  ".RData"
# )
# 
# load(file_in)

f_load_if_absent_csv(dir_in = dir_in,
                     name_env_variable = "coldata",
                     match_str = "METADATA")

cts <- f_read_csv(dir = dir_in,
                  csv_file = f_list_files(dir = dir_in,
                                          match_str = "P0_CTS_gen",
                                          file_type = "csv"
                                          )
                  )

cts_filtered <- f_read_csv(dir = dir_in,
                           csv_file = f_list_files(dir = dir_in,
                                                   match_str = "P1B1_cts_filtered",
                                                   file_type = "csv"
                           )
)

# Calculate rowsum before/after filtering
cts_sum <- setNames(as.data.frame(colSums(cts)), "cts_filterd_sum")
cts_filtered_sum <- setNames(as.data.frame(colSums(cts_filtered)), "cts_filterd_sum")

# Add columns to coldata
coldata <- cbind(coldata, cts_sum, cts_filtered_sum) 

####################### Plot info 
###################### 

colsum_df <- data.frame(sample = as.factor(rownames(coldata)),
                        reads_Pe = sum_Pe, 
                        reads_Sp = sum_Sp, 
                        reads_all = coldata_sorted$reads
                        )


# Safe for later usage
write.csv(colsum_df, paste0(dir_out, "read_info.csv"), row.names = FALSE)

#coldata_reads <- cbind(coldata_sorted, colsum_df)

# 4.3 Convert the dataframe from wide to long format
colsum_long_all <- tidyr::pivot_longer(colsum_df, 
                               cols = starts_with("reads_"), 
                               names_to = "source", 
                               values_to = "read_count")

rm(colsum_df)

# Filter rows without 'reads_unm'
colsum_long_unm <- colsum_long_all %>% 
  filter(!(source %in% "reads_all"))

colsum_long_org <- colsum_long_unm %>% 
  filter(!(source %in% "reads_unm"))

colsum_long_reads <- colsum_long_all %>% 
  filter(source %in% c("reads_all"))

data_plot <- list(colsum_long_unm, colsum_long_org, colsum_long_reads)
plot_names <- c("bargraph_unmapped", "bargraph_mapped", "bargraph_reads")

# 5. Draw and save staged bar graphs ----
# 5.1 Function for boxplot
f_bar_reads <- function(df){
g <- ggplot(df, aes(
  x = factor(sample, levels = sample_order), 
  y = read_count, 
  fill = source)) +
  geom_col()+
  scale_fill_manual(values = c("#E69F00", "#009E73", "#5C5C5C")) +
  # define the theme of the boxplot
  theme_bw() +  # make the bg white
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_blank(), # remove background, frame
        axis.line = element_line(colour = "black"))+
  # label the axises 
  xlab("sample") +                
  ylab("mapped reads") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

return(g)
}

# 5.2 Apply the function over the list of data frames
plots <- lapply(data_plot, f_bar_reads)

# 5.3 Combine plots and their names into a named list
plot_list <- setNames(plots, plot_names)

# 5.4 Save each plot with its respective name
lapply(names(plot_list), function(plot_name) {
  plot <- plot_list[[plot_name]]
  ggsave(paste0(dir_out, "reads_", plot_name, ".svg"), plot, 
         device = "svg",
         width = 16,
         height = 5)
})

