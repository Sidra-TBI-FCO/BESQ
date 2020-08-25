
# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA")                                                                   # Setwd to location were output files have to be saved.

source("../tools/ipak.function.R")
ipak(c("ggplot2", "ggpubr", "dplyr"))

# Set parameters


# Load data
load("../tools/Selected.pathways.3.4.RData")
load(paste0("./Data/expression matrix/BRCA_gene_RNAseq_normalized_Tissue_Type_and_whitelists_filtered.Rdata"))
load("./Analysis/Sample_annotations.Rdata")
load("./Analysis/t_test_between_ethnicities/C23_significant_pathways_BasalMyo_All.Rdata")
patient_data = read.csv("./Data/Clinical Data/patient.csv", stringsAsFactors = FALSE)
patient_data = patient_data[-c(1, 2),]
#load(paste0("./Analysis/t_test_between_ethnicities/C14_Diff_", Gene_set, "_between_ICR_IMS_Ethnicity.Rdata"))

# Set parameters
Ethnicity_variable = "Assigned_Ethnicity_simplified" #"Ethnicity_reported_simple"
variables = significant_pathways
#variables = rownames(ES)

# Analysis
if(Ethnicity_variable == "Assigned_Ethnicity_simplified"){
  Clinical_data_ann = Clinical_data_ann[-which(Clinical_data_ann[,Ethnicity_variable] %in% c("Unclear", "Asian")),]
}

#Clinical_data_ann = Clinical_data_ann[which(Clinical_data_ann$IMS_Mathews == "BasalMyo"),]
table(Clinical_data_ann$IMS_Mathews, exclude = NULL)
Clinical_data_ann$IMS_Mathews = as.character(Clinical_data_ann$IMS_Mathews)
Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann$IMS_Mathews)),]

Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann[, Ethnicity_variable])),]

Clinical_data_ann$HML.ICR.Cluster = factor(Clinical_data_ann$HML.ICR.Cluster, levels = c("ICR High", "ICR Medium", "ICR Low"))

Clinical_data_ann$Estrogen_receptor = patient_data$er_status_by_ihc[match(Clinical_data_ann$bcr_patient_barcode,
                                                                          patient_data$bcr_patient_barcode)]

table(Clinical_data_ann$Estrogen_receptor, exclude = NULL)
Clinical_data_ann$Estrogen_receptor[which(Clinical_data_ann$Estrogen_receptor %in% c("[Not Evaluated]", "Indeterminate"))] = NA
Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann$Estrogen_receptor)),]
Clinical_data_ann$Category = NA
Clinical_data_ann$Category[which(Clinical_data_ann$Estrogen_receptor == "Positive")] = "ER positive"
Clinical_data_ann$Category[which(Clinical_data_ann$IMS_Mathews == "BasalMyo")] = "BasalMyo"

genes = Selected.pathways$`[IPA] ErbB Signaling`
genes = genes[which(genes %in% rownames(filtered.norm.RNAseqData))]
Expression_data = filtered.norm.RNAseqData[genes,]
Expression_data = log(Expression_data +1, 2)
dim(Expression_data)

results_df = data.frame(Gene = rownames(Expression_data), 
                        p.value = 0, p.value.fdr = 0, fc = 0, 
                        mean_BasalMyo = 0, mean_ER_positive = 0,
                        CI_lower = 0, CI_upper = 0)
data = data.frame(Patient_ID = substring(colnames(Expression_data), 1, 12),
                  Category = NA,
                  Expression = NA)
data$Category = Clinical_data_ann$Category[match(data$Patient_ID, Clinical_data_ann$bcr_patient_barcode)]

data = data[-which(is.na(data$Category)),]

i=1
for (i in 1:nrow(Expression_data)){
  gene = rownames(Expression_data)[i]
  data$Expression = Expression_data[gene,][match(data$Patient_ID, substring(colnames(Expression_data), 1, 12))]
  t_test = t.test(data$Expression[which(data$Category == "BasalMyo")], data$Expression[which(data$Category == "ER positive")])
  results_df$mean_BasalMyo[which(results_df$Gene == gene)] = t_test$estimate[1]
  results_df$mean_ER_positive[which(results_df$Gene == gene)] = t_test$estimate[2]
  results_df$p.value[which(results_df$Gene == gene)] = t_test$p.value
  results_df$CI_lower[which(results_df$Gene == gene)] = t_test$conf.int[1]
  results_df$CI_upper[which(results_df$Gene == gene)] = t_test$conf.int[2]
  results_df$fc[which(results_df$Gene == gene)] = t_test$estimate[1] -  t_test$estimate[2]
}
results_df$p.value.fdr = p.value.fdr <- p.adjust(p = results_df$p.value,method = "fdr",n = nrow(results_df))

dir.create("./Analysis/Reviewer_comment_Erbb", showWarnings = FALSE)
write.table(results_df, file = "./Analysis/Reviewer_comment_Erbb/results_ErbB_signaling.tsv", sep = "\t", quote = FALSE,
            col.names = NA)
write.csv(results_df, file = "./Analysis/Reviewer_comment_Erbb/results_ErbB_signaling.csv")
