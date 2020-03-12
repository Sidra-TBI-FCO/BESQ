
# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA")                                                                   # Setwd to location were output files have to be saved.

source("../tools/ipak.function.R")

required.packages = c("ggplot2")
ipak(required.packages)

# Load data
load("./Data/expression matrix/BRCA_gene_RNAseq_normalized_Tissue_Type_and_whitelists_filtered.Rdata")
RNASeq.QN.LOG2 = log(filtered.norm.RNAseqData +1, 2)
load("./Analysis/Sample_annotations.Rdata")
load("Analysis/ssGSEA_scores/Selected_pathways_ES_scores.Rdata")

# Set parameters
Ethnicity = "Black"

# Analysis
ES = as.data.frame(t(ES))
Clinical_data_ann = Clinical_data_ann[which(Clinical_data_ann$IMS_Mathews == "BasalMyo"),]
Clinical_data_ann = Clinical_data_ann[which(Clinical_data_ann$Assigned_Ethnicity_simplified %in% Ethnicity),]

plot_df = data.frame(Sample = Clinical_data_ann$bcr_patient_barcode, ICR_ES = Clinical_data_ann$ICR_ES)
plot_df$AMPK_sig = ES$`[IPA] AMPK Signaling`[match(plot_df$Sample, substring(rownames(ES), 1, 12))]

plot = ggplot(plot_df, aes(x = AMPK_sig, y = ICR_ES)) +
  geom_point() + theme_bw()

dev.new()
plot(plot)
