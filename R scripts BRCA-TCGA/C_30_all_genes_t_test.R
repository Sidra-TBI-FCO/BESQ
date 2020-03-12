
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

# Analysis
Clinical_data_ann = Clinical_data_ann[which(Clinical_data_ann$Assigned_Ethnicity_simplified %in% c("Black", "White")),]
#Clinical_data_ann = Clinical_data_ann[which(Clinical_data_ann$HML.ICR.Cluster %in% c("ICR Low")),]
Clinical_data_ann = Clinical_data_ann[which(Clinical_data_ann$IMS_Mathews %in% c("BasalMyo")),]

colnames(RNASeq.QN.LOG2) = substring(colnames(RNASeq.QN.LOG2), 1, 12)
RNASeq.QN.LOG2 = RNASeq.QN.LOG2[,Clinical_data_ann$bcr_patient_barcode]

Genes = rownames(RNASeq.QN.LOG2)
N.Genes = length(Genes)

results = data.frame(Gene = Genes, p_value = NA, p_value_FDR = NA, mean_W = NA, mean_B = NA)

i = 2
for (i in 1:N.Genes){
  gene = Genes[i]
  df_test = data.frame(Patient_ID = colnames(RNASeq.QN.LOG2), 
                       Expression = RNASeq.QN.LOG2[which(rownames(RNASeq.QN.LOG2) == gene),], 
                       Ethnicity = NA)
  df_test$Ethnicity = Clinical_data_ann$Assigned_Ethnicity_simplified[match(df_test$Patient_ID,
                                                                            Clinical_data_ann$bcr_patient_barcode)]
  t_test = t.test(df_test$Expression[which(df_test$Ethnicity == "White")],
                  df_test$Expression[which(df_test$Ethnicity == "Black")],
                  paired = FALSE,
                  conf.level = 0.95)
  results$mean_W[which(results$Gene == gene)] = mean(df_test$Expression[which(df_test$Ethnicity == "White")])
  results$mean_B[which(results$Gene == gene)] = mean(df_test$Expression[which(df_test$Ethnicity == "Black")])
  results$p_value[which(results$Gene == gene)] = t_test$p.value
  results$p_value_FDR[which(results$Gene == gene)] = p.adjust(t_test$p.value, method = "fdr", n = nrow(results))
}

results$log_fold_change = results$mean_B - results$mean_W

dir.create("./Analysis/C_30_all_genes_t_test", showWarnings = FALSE)
save(results, file = "./Analysis/C_30_all_genes_t_test/diff_expressed_genes_between_B_W_in_All_ICR_groups.Rdata")
write.csv(results, file = "./Analysis/C_30_all_genes_t_test/diff_expressed_genes_between_B_W_in_All_ICR_groups.csv")

results_sig_0.05 = results[which(results$p_value < 0.05),]
write.csv(results_sig_0.05, file = "./Analysis/C_30_all_genes_t_test/sig_diff_0.05_genes_between_B_W_in_All_ICR_groups.csv")
