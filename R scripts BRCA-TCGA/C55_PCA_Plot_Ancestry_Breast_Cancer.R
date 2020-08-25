
# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA")                                                                   # Setwd to location were output files have to be saved.
source("../tools/ipak.function.R")
ipak(c("ggplot2", "ggpubr"))

# Load data
load("./Analysis/Sample_annotations.Rdata")
UCSF_Ancestry_Calls = read.csv("./Data/Ancestry data/UCSF_Ancestry_Calls.csv", stringsAsFactors = FALSE)

# Plot
plot_df = Clinical_data_ann[, c("bcr_patient_barcode","Assigned_Ethnicity_simplified")]
plot_df$PC1 = UCSF_Ancestry_Calls$PC1[match(plot_df$bcr_patient_barcode, UCSF_Ancestry_Calls$Patient_ID)]
plot_df$PC2 = UCSF_Ancestry_Calls$PC2[match(plot_df$bcr_patient_barcode, UCSF_Ancestry_Calls$Patient_ID)]

plot_df = plot_df[-which(is.na(plot_df$PC1)),]
plot_df = plot_df[-which(plot_df$Assigned_Ethnicity_simplified == "Unclear"),]

plot = ggplot(plot_df, aes(x = PC1, y = PC2, color = Assigned_Ethnicity_simplified)) +
  geom_point(size = 0.4) +
  scale_color_manual(values = c("White" = "#FF7E00", 
                                "Black" = "#066C3C", 
                                "Asian" = "#7F2D80")) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.x = element_text(color = "black", size = 12),
        axis.title.y = element_text(color = "black", size = 12),
        aspect.ratio = 1/1,
        legend.position = "none")

dir.create("./Figures/C55_PCA_Plot_Ancestry_Breast", showWarnings = FALSE)
png("./Figures/C55_PCA_Plot_Ancestry_Breast/PCA_Plot_Ancestry_breast.png",
    res = 600, width = 3, height = 3.5, units = "in")
plot(plot)
dev.off()
