
#Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA")                                                                   # Setwd to location were output files have to be saved.

source("../tools/ipak.function.R")

ipak(c("ggplot2", "ggpubr", "dplyr"))

# load data
load("./Analysis/Sample_annotations.Rdata")
load("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/3_DataProcessing/External/BinaryMatrix.From.Michele.RData")

# Set parameters
Ethnicity_variable = "Assigned_Ethnicity_simplified"

# Analysis
Clinical_data_ann = Clinical_data_ann[-which(Clinical_data_ann[, Ethnicity_variable] %in% c("Other", "Asian", "Hispanic", "Unclear")),]
Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann$HML.ICR.Cluster)),]
Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann[, Ethnicity_variable])),]
Clinical_data_ann = Clinical_data_ann[which(Clinical_data_ann$IMS_Mathews == "BasalMyo"),]
Clinical_data_ann$Assigned_Ethnicity_simplified = factor(Clinical_data_ann$Assigned_Ethnicity_simplified, levels = c("White", "Black"))

binary_matrix = binary_matrix[which(substring(rownames(binary_matrix), 1, 12) %in% Clinical_data_ann$bcr_patient_barcode),1:470]
Clinical_data_ann = Clinical_data_ann[which(Clinical_data_ann$bcr_patient_barcode %in% substring(rownames(binary_matrix), 1, 12)),]

binary_matrix_sub = binary_matrix[, which(colSums(binary_matrix) > 2)]

results = data.frame(Gene = colnames(binary_matrix_sub),
                     N_mut_white = NA,
                     N_WT_white = NA,
                     N_mut_black = NA,
                     N_WT_black = NA,
                     p_val_chisq = NA)

i=2
for (i in 1:ncol(binary_matrix_sub)){
  gene = colnames(binary_matrix_sub)[i]
  df = data.frame(Patient = Clinical_data_ann$bcr_patient_barcode,
                  Ethnicity = Clinical_data_ann$Assigned_Ethnicity_simplified,
                  Mutation = NA)
  df$Mutation = binary_matrix_sub[, gene][match(df$Patient, substring(rownames(binary_matrix_sub), 1, 12))]
  tbl = table(df$Ethnicity, df$Mutation)
  chisq = chisq.test(tbl)
  results$p_val_chisq[which(results$Gene == gene)] = chisq$p.value
  results$N_mut_black[which(results$Gene == gene)] = tbl["Black", "1"]
  results$N_WT_black[which(results$Gene == gene)] = tbl["Black", "0"]
  results$N_mut_white[which(results$Gene == gene)] = tbl["White", "1"]
  results$N_WT_white[which(results$Gene == gene)] = tbl["White", "0"]
}

results$percent_mutated_in_white = (results$N_mut_white/128)*100
results$percent_mutated_in_black = (results$N_mut_black/58)*100
