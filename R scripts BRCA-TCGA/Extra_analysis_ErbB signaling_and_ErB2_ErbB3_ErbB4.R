

# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA")                                                                   # Setwd to location were output files have to be saved.

source("../tools/ipak.function.R")

ipak(c("ggplot2", "ggpubr", "dplyr"))

# Set parameters
Gene_set = "Selected.pathways" #"Bindea_ORIG"  #"Selected_pathways"
Stages = ""
Ethnicity_variable = "Assigned_Ethnicity_simplified" #"Ethnicity_reported_simple"

# Load data
load(paste0("./Analysis/ssGSEA_scores/v2_", Gene_set, "_ES_scores.Rdata"))
load("./Analysis/Sample_annotations.Rdata")

if(Ethnicity_variable == "Assigned_Ethnicity_simplified"){
  Clinical_data_ann = Clinical_data_ann[-which(Clinical_data_ann[,Ethnicity_variable] %in% c("Unclear", "Asian")),]
}

Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann[, Ethnicity_variable])),]
Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann$IMS_Mathews)),]
Clinical_data_ann = Clinical_data_ann[which(Clinical_data_ann$IMS_Mathews == "BasalMyo"),]
Clinical_data_ann = Clinical_data_ann[which(Clinical_data_ann$Assigned_Ethnicity_simplified == "Black"),]

# Scatterplots
df_plot = Clinical_data_ann[, c("bcr_patient_barcode", "IMS_Mathews", Ethnicity_variable)]
df_plot$ErbB_signaling = ES["[IPA] ErbB Signaling",][match(df_plot$bcr_patient_barcode, substring(colnames(ES), 1, 12))]
df_plot$ErbB2_ErbB3_ErbB4 = ES["ERBB2, ERBB3, and ERBB4",][match(df_plot$bcr_patient_barcode, substring(colnames(ES), 1, 12))]


plot = ggplot(df_plot, aes(x = ErbB_signaling, y = ErbB2_ErbB3_ErbB4)) +
  geom_point(size = 0.4) +
  stat_cor(method = "pearson", size = 5, label.y = -0.2, label.x = 0.2) +
  geom_smooth(method="lm") +
  xlab("ErbB signaling") +
  ylab("ErbB2+ErbB3+ErbB4") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 17, colour = "black"),
        axis.title.y = element_text(size = 17, colour = "black"),
        axis.text.x = element_text(size = 17, colour = "black"),
        axis.text.y = element_text(size = 17, colour = "black"))

dev.new()  
plot(plot)
