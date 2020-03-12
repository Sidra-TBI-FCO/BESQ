# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA") 

source("~/Dropbox (TBI-Lab)/DB-LAB/bioinformatics tools/R scripts/ipak.function.R")

required.packages <- c("corrplot", "stringr", "ggplot2", "ggpubr")
ipak(required.packages)

# Set parameters
variable_of_interest = "StromalScore" #"ImmuneScore" #"StromalScore" #"ESTIMATEScore"

# Create directories
dir.create("./Figures/C_20_ESTIMATE", showWarnings = FALSE)

# Load data
load("./Analysis/ESTIMATE/TCGA_BRCA_ESTIMATE_scores.Rdata")
load("./Analysis/Sample_annotations.Rdata")

ESTIMATE = as.data.frame(ESTIMATE)

ESTIMATE$Ethnicity = Clinical_data_ann$Assigned_Ethnicity_simplified[match(substring(rownames(ESTIMATE), 1, 12), Clinical_data_ann$bcr_patient_barcode)]
ESTIMATE$IMS_Mathews = Clinical_data_ann$IMS_Mathews[match(substring(rownames(ESTIMATE), 1, 12), Clinical_data_ann$bcr_patient_barcode)]

df_plot = ESTIMATE[, c("Ethnicity", "IMS_Mathews", variable_of_interest)]
colnames(df_plot) = c("Ethnicity", "IMS_Mathews", "variable_of_interest")

df_plot = df_plot[which(df_plot$Ethnicity %in% c("White", "Black")),]

my_comparisons = list(c("Black", "White"))
plot = ggplot(df_plot, aes(x = Ethnicity, y = variable_of_interest)) +
  geom_boxplot(outlier.shape = NA) +
  ylab(variable_of_interest) +
  geom_jitter(size = 0.2, width = 0.2, aes(color = Ethnicity)) +
  scale_color_manual(values = c("#FF7E00", "#066C3C")) +
  theme_bw() +
  facet_grid(.~IMS_Mathews) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test")

png(paste0("./Figures/C_20_ESTIMATE/", variable_of_interest, "_White_vs_Black.png"), res = 600, height = 4, width = 8, units = "in")
plot(plot)
dev.off()
