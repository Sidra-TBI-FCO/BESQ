
# Set-up environment 
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA")
source("../tools/ipak.function.R")
required.packages <- c("ggplot2")
ipak(required.packages)  

# Load data
load("./Analysis/Sample_annotations.Rdata")
load("./Analysis/ssGSEA_scores/Selected_pathways_ES_scores.Rdata")
#

# Set parameters
IMS = "Non_BasalMyo"
Pathway = "[IPA] AMPK Signaling" #"[HM] TGF beta signaling" #"[IPA] ErbB2 ErbB3 Signaling" #"[HM] Glycolysis" #"[IPA] AMPK Signaling"
ICR = "All"
Ethnicity = "White"

load(paste0("./Analysis/C_26_Selected_pathway_classification/", IMS,"_median_[IPA] AMPK Signaling_Blackclassification.Rdata"))
# Classification
#Clinical_data_ann = Clinical_data_ann[which(Clinical_data_ann$IMS_Mathews == IMS),]
Clinical_data_ann = Clinical_data_ann[which(Clinical_data_ann$Assigned_Ethnicity_simplified %in% Ethnicity),]
if(IMS == "Non_BasalMyo"){
  Clinical_data_ann = Clinical_data_ann[-which(Clinical_data_ann$IMS_Mathews == "BasalMyo"),]
}
#Clinical_data_ann = Clinical_data_ann[which(Clinical_data_ann$HML.ICR.Cluster %in% ICR),]

ES = ES[,which(substring(colnames(ES), 1, 12) %in% Clinical_data_ann$bcr_patient_barcode)]

df = data.frame(Sample = colnames(ES),
                Variable = ES[Pathway,],
                Category = NA)

if(Ethnicity == "Black"){
  med = median(df$Variable)
}else{med=med}
df$Category[which(df$Variable < med)] = "Pathway Low"
df$Category[which(df$Variable >= med)] = "Pathway High"

# Save data
dir.create("./Analysis/C_26_Selected_pathway_classification", showWarnings = FALSE)
save(df, file = paste0("./Analysis/C_26_Selected_pathway_classification/", IMS, "_", Pathway, "_", Ethnicity,"_binary_classification.Rdata"))
save(med, file = paste0("./Analysis/C_26_Selected_pathway_classification/", IMS, "_median_", Pathway,"_", Ethnicity, "classification.Rdata"))
#df$Category = NA
#first = unname(quantile(df$Variable, probs = seq(0, 1, length.out = 4))[2])
#second = unname(quantile(df$Variable, probs = seq(0, 1, length.out = 4))[3])
#df$Category[which(df$Variable < first)] = "Pathway Low"
#df$Category[which(df$Variable >= first & df$Variable < second)] = "Pathway Medium"
#df$Category[which(df$Variable >= second)] = "Pathway High"

#save(df, file = paste0("./Analysis/C_26_Selected_pathway_classification/BasalMyo_", Pathway, "_", Ethnicity, "_tertiles_classification.Rdata"))

