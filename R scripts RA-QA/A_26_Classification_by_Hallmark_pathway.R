
# Set-up environment 
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-RA-QA")
source("../tools/ipak.function.R")
required.packages <- c("ggplot2")
ipak(required.packages)  

# Load data
load("./Analysis/Sample_annotations.Rdata")
load("./Analysis/ssGSEA_scores/Selected_pathways_ES_scores.Rdata")

# Set parameters
Pathway = "[IPA] AMPK Signaling" #"[HM] TGF beta signaling" #"[IPA] ErbB2 ErbB3 Signaling" #"[HM] Glycolysis" #"[IPA] AMPK Signaling"
ICR = "All"
Ethnicity = "Arab"

# Classification
#Clinical_data_ann = Clinical_data_ann[which(Clinical_data_ann$IMS_Mathews == IMS),]
Clinical_data_ann = Clinical_data_ann[which(Clinical_data_ann$Ethnicity %in% Ethnicity),]
Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann$HL.ICR.Cluster)),]

ES = ES[,which(colnames(ES) %in% Clinical_data_ann$Sample_ID)]

df = data.frame(Sample = colnames(ES),
                Variable = ES[Pathway,],
                Category = NA)


med = median(df$Variable)

df$Category[which(df$Variable < med)] = "Pathway Low"
df$Category[which(df$Variable >= med)] = "Pathway High"

# Save data
dir.create("./Analysis/A_26_Selected_pathway_classification", showWarnings = FALSE)
save(df, file = paste0("./Analysis/A_26_Selected_pathway_classification/", Pathway, "_", Ethnicity,"_binary_classification.Rdata"))
#save(med, file = paste0("./Analysis/C_26_Selected_pathway_classification/BasalMyo_median_", Pathway,"_", Ethnicity, "classification.Rdata"))


