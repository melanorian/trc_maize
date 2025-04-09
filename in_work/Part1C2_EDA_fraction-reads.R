# This script allows the selection of specific samples from the data set
# The output is a .csv file including the raw counts in the samples of interest

# clean environment 
rm(list = ls())

while (dev.cur() != 1) {
  dev.off()}

# 0. base set up install/load libraries, setwd, read files ----
# 0.1 load libraries
library(tidyr)
library(ggplot2)

# 1.1 load def_variables.R and def_functions.R
source("/media/wolf2/TKI2019/melanie/Code-transcriptomics-analysis/Pipeline/def_variables.R")
source("/media/wolf2/TKI2019/melanie/Code-transcriptomics-analysis/Pipeline/def_functions.R")

# 1.2 Set wd
setwd(dir_base)

# 1.3 in direcory
dir_in <- paste0(dir_base, "DDS/")

# 1.4 out directory
dir_eda <- paste0(dir_base, "EDA_basic/")

name_out <- f_make_name(select_genes = select_genes,
                        remove_samples = remove_samples)

dir_out <- paste0(dir_eda, name_out, "/")

if (!file.exists(dir_eda)) {
  dir.create(dir_eda)
}

if (!file.exists(dir_out)) {
  dir.create(dir_out)
}

# Load info on reads
df <- f_read_csv(dir = dir_out, csv_file = "read_info.csv")

# 6. Draw scatter plot ----
# 4. Get the total depth for each sample
graph <- "disease-progression"

df$reads_map <- df$reads_Pe + df$reads_Sp

# Calculate the percentage of reads mapping to Pe and Sp for each sample (total reads or mapped reads)
df$pe_percentages <- (df$reads_Pe / df$reads_map) * 100
df$sp_percentages <- (df$reads_Sp / df$reads_map) * 100

# simpler df
keep_cols <- names(df)[grep(paste0("pe_percentages|sp_percentages", "|", factor_experiment), names(df), invert = F)]
df_simpl <- df[, keep_cols]

# Assuming percentage_results is your data frame containing Pe and Sp percentages
df_long <- df_simpl %>%
  tidyr::pivot_longer(cols = c(pe_percentages, sp_percentages), names_to = "Organism", values_to = "Percentage")

# Convert dpi to a factor with the desired levels
df_long[[factor_experiment]] <- as.factor(df_long[[factor_experiment]])

# 1.3 Grouped jitterplot with lines 
# This requires additional statistics for plotting (sd, )
df_sum <- df_long %>%
  mutate(Organism = factor(Organism)) %>% 
  group_by(!!sym(factor_experiment), Organism) %>%
  summarise(
    sd = sd(Percentage),
    Percentage = mean(Percentage)
  )

g <- ggplot(df_long, aes(x = !!sym(factor_experiment) , y= Percentage, color = Organism)) +
  geom_jitter(position = position_jitter(0.2), size = 5, alpha = 0.5) + 
  geom_line(aes(group = Organism), data = df_sum, linejoin = "mitre", linemitre = 10, lineend = "round") +
  geom_errorbar(aes(ymin = Percentage-sd, ymax = Percentage+sd), data = df_sum, width = 0.2)+
  scale_color_manual(values = c("#E69F00", "#009E73")) +
  geom_point(data = df_sum, size = 7, alpha = 0.9) +
  theme(legend.position = "top") +
  theme_classic() +
  ylab("Percentage of total reads mapping to Organism")+
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 20))



# add stats to g

filter <- df_long[grep("pe_", df_long$Organism, invert = F), ] # filter df

# ANOVA to analyse whether there are differences between the groups
# H0 All means are the same, H1 at least one mean is different
formula <- as.formula(paste("Percentage ~", factor_experiment))
anova_reads <- aov(formula, data = filter)

sum_reads_aov <- summary(anova_reads)

# TUKEY - which groups are significantly different 
TUKEY <- TukeyHSD(anova_reads)

# Summarize TUKEY test results for saving & plotting

# 4.3.1 write a funktion to generate labels for the test results
generate_label_df <- function(TUKEY, variable){
  
  # Extract labels and factor levels from Tukey post-hoc results
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  # order labels according to boxplot:
  Tukey.labels$dpi = rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$dpi) , ]
  return(Tukey.labels)
}

# 4.3.2 Apply the function on the df
LABELS <- generate_label_df(TUKEY , paste0(factor_experiment))
names(LABELS) <- c("Letters",paste0(factor_experiment))

# Merge the label data with the original dataframe for plotting
df_long_labels <- merge(df_long, LABELS, 
                             by.x = as.character(factor_experiment), 
                             by.y = "row.names", 
                             all.x = TRUE
                             )

df_long_labels$yloc <- 105
  
# Replace Letters with empty string for rows where Organism starts with "sp"
df_long_labels[startsWith(df_long_labels$Organism, "sp"), "Letters"] <- ""

# draw boxplot
g2 <- ggplot(df_long_labels, aes(x = !!sym(factor_experiment) , y= Percentage, color = Organism)) +
  geom_jitter(position = position_jitter(0.2), size = 5, alpha = 0.5) + 
  geom_line(aes(group = Organism), data = df_sum, linejoin = "mitre", linemitre = 10, lineend = "round") +
  geom_errorbar(aes(ymin = Percentage-sd, ymax = Percentage+sd), data = df_sum, width = 0.2)+
  scale_color_manual(values = c("#E69F00", "#009E73")) +
  geom_point(data = df_sum, size = 7, alpha = 0.9) +
  theme(legend.position = "top") +
  theme_classic() +
  geom_text(data = df_long_labels, aes(y = yloc, label=Letters)) +
  ylab("Percentage of total reads mapping to Organism")+
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        plot.margin = margin(50, 10, 10, 10, "cm"))

print(g2)

ggsave(paste0(dir_out, "read-fractions",".svg"), g, 
       device = "svg",
       width = 13,
       height = 8)

dev.off()
