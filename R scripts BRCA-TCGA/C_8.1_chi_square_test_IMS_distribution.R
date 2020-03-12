

# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA")
source("../tools/ipak.function.R")
ipak(c("dplyr", "ggplot2"))

#load data
load("../RNAseq-RA-QA/Analysis/Sample_annotations.Rdata")
Clinical_data_ann_RA_QA = Clinical_data_ann

load("./Analysis/Sample_annotations.Rdata")

Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann$IMS_Mathews)),]
Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann$Assigned_Ethnicity_simplified)),]
Clinical_data_ann$IMS_Mathews = as.character(Clinical_data_ann$IMS_Mathews)
Clinical_data_ann$Assigned_Ethnicity_simplified = as.character(Clinical_data_ann$Assigned_Ethnicity_simplified)

colnames(Clinical_data_ann_RA_QA)[which(colnames(Clinical_data_ann_RA_QA) == "Ethnicity")] = "Assigned_Ethnicity_simplified"
Clinical_data_ann = rbind(Clinical_data_ann[, c("IMS_Mathews", "Assigned_Ethnicity_simplified")], Clinical_data_ann_RA_QA[,c("IMS_Mathews","Assigned_Ethnicity_simplified")])

Clinical_data_ann = Clinical_data_ann[-which(Clinical_data_ann$Assigned_Ethnicity_simplified %in% c("Unclear", "Black", "Asian", "Other")),]

# Analysis
subtypes = unique(Clinical_data_ann$IMS_Mathews)
results = data.frame(Subtype= subtypes, X_squared = NA, p_value = NA)

i = 5
for (i in 1:length(subtypes)){
  subtype = subtypes[i]
  analysis = Clinical_data_ann
  analysis$IMS_Mathews[-which(analysis$IMS_Mathews == subtype)] = "Other IMS"
  analysis$IMS_Mathews[which(analysis$IMS_Mathews == subtype)] = "Interest IMS"
  tbl = table(analysis$Assigned_Ethnicity_simplified, analysis$IMS_Mathews)
  result = chisq.test(tbl)
  results$X_squared[which(results$Subtype == subtype)] = result$statistic
  results$p_value[which(results$Subtype == subtype)] = result$p.value
}
ethnicities = paste(unique(analysis$Assigned_Ethnicity_simplified), collapse = "_")

dir.create("./Analysis/C8.1_Chi_squared_test_results", showWarnings = FALSE)
write.csv(results, file = paste0("./Analysis/C8.1_Chi_squared_test_results/", ethnicities, "_IMS_Mathews_Chi_squared_result.csv"))

