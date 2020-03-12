#
# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-RA-QA") 

source("~/Dropbox (TBI-Lab)/DB-LAB/bioinformatics tools/R scripts/ipak.function.R")

required.packages <- c("corrplot", "stringr", "ggplot2", "ggpubr")
ipak(required.packages)

# Set parameters
variable_of_interest = "ESTIMATEScore" #"ImmuneScore" #"StromalScore" #"ESTIMATEScore"

# Create directories
dir.create("./Figures/A_20_ESTIMATE", showWarnings = FALSE)

# Load data
load("./Analysis/ESTIMATE/TCGA_BRCA_ESTIMATE_scores.Rdata")
load("./Analysis/Sample_annotations.Rdata")

ESTIMATE = as.data.frame(ESTIMATE)

ESTIMATE$Ethnicity = Clinical_data_ann$Ethnicity[match(rownames(ESTIMATE), Clinical_data_ann$Sample_ID)]
ESTIMATE$IMS_Mathews = Clinical_data_ann$IMS_Mathews[match(rownames(ESTIMATE), Clinical_data_ann$Sample_ID)]

df_plot = ESTIMATE[, c("Ethnicity", "IMS_Mathews", variable_of_interest)]
colnames(df_plot) = c("Ethnicity", "IMS_Mathews", "variable_of_interest")

plot = ggplot(df_plot, aes(x = Ethnicity, y = variable_of_interest)) +
  geom_boxplot(outlier.shape = NA) +
  ylab(variable_of_interest) +
  geom_jitter(size = 0.2, width = 0.2, aes(color = Ethnicity)) +
  scale_color_manual(values = c("violet", "#7F2D80", "lightblue")) +
  theme_bw() +
  facet_grid(.~IMS_Mathews)

png(paste0("./Figures/A_20_ESTIMATE/", variable_of_interest, ".png"), res = 600, height = 4, width = 11, units = "in")
plot(plot)
dev.off()
