
# Set up environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA")                                                                   # Setwd to location were output files have to be saved.

source("../tools/ipak.function.R")

# Set parameters
Ethnicity_variable = "Assigned_Ethnicity_simplified" #"Ethnicity_reported_simple"
IMS_group = c("BasalMyo")
ICR_group = "All"

# Load data
load("./Analysis/Sample_annotations.Rdata")
load("./Analysis/ssGSEA_scores/Selected_pathways_ES_scores.Rdata")
Survival_gx_boost = read.csv("./Analysis/Gxboost_trees_analysis_Raghvendra/Results/Feature_Importance_For_Survival.csv", stringsAsFactors = FALSE)
Survival_gx_boost$Features = paste("[", Survival_gx_boost$Features, sep = "")
Survival_gx_boost$Features = gsub("\\:", "] ", Survival_gx_boost$Features)
Survival_gx_boost$Features = gsub("\\_", " ", Survival_gx_boost$Features)

# Analysis
Clinical_data_ann = Clinical_data_ann[which(Clinical_data_ann$IMS_Mathews == IMS_group),]
Clinical_data_ann = Clinical_data_ann[which(Clinical_data_ann$Assigned_Ethnicity_simplified %in%
                                              c("White", "Black")),]
if(ICR_group == "All"){}else{
  Clinical_data_ann = Clinical_data_ann[which(Clinical_data_ann$HML.ICR.Cluster == ICR_group),]
}
colnames(ES) = substring(colnames(ES), 1, 12)

ES = ES[,which(colnames(ES) %in% Clinical_data_ann$bcr_patient_barcode)]
ES = t(ES)

results_df = data.frame(Pathway = colnames(ES), Mean_white = NA, Mean_black = NA,
                        p_value = NA, CI_lower = NA, CI_upper = NA)
results_df$Pathway = as.character(results_df$Pathway)

i = 1
for (i in 1:nrow(results_df)){
  pathway = results_df$Pathway[i]
  df = data.frame(Sample = rownames(ES), Variable = ES[,which(colnames(ES) == pathway)], Ethnicity = NA)
  df$Ethnicity = Clinical_data_ann$Assigned_Ethnicity_simplified[match(df$Sample, Clinical_data_ann$bcr_patient_barcode)]
  result = t.test(Variable ~ Ethnicity, df)
  results_df[i, 2:6] = c(result$estimate["mean in group White"], result$estimate["mean in group Black"],
                         result$p.value, result$conf.int[1], result$conf.int[2])
}

results_df$FDR = p.adjust(p = results_df$p_value, method = "fdr",n = nrow(results_df))

significant_pathways = results_df$Pathway[which(results_df$p_value < 0.05)]
Survival_gx_boost = Survival_gx_boost[which(Survival_gx_boost$Features %in% significant_pathways),]
significant_pathways = Survival_gx_boost$Features[1:8]

save(significant_pathways, file = paste0("./Analysis/t_test_between_ethnicities/C23_significant_pathways_and_importance_gxboost_tree_survival_",
                                         IMS_group, "_", ICR_group, ".Rdata"))
