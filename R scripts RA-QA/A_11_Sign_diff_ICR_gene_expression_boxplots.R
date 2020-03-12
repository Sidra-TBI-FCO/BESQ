

# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-RA-QA")                                                                   # Setwd to location were output files have to be saved.

source("../tools/ipak.function.R")

ipak(c("ggplot2", "ggpubr"))

# Set parameters

# Load data
load("./Analysis/Sample_annotations.Rdata")
load("./Data/expression matrix/CTA_TNBC_BRCA_Qatar_RNASeq.tumor.tissue.filtered.EDAseq.QN.Rdata")
load("../tools/ICR_genes.RData")

dir.create("./Figures/C11_Boxplot_ICR_gene_across_ethnicities", showWarnings = FALSE)

# Analysis

# Delete NA
Clinical_data_ann = Clinical_data_ann[-which(Clinical_data_ann$Ethnicity == "Other"),]
Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann$IMS)),]

ICR_clusters = names(table(Clinical_data_ann$HL.ICR.Cluster))
IMS_categories = names(table(Clinical_data_ann$IMS))
Ethicities = names(table(Clinical_data_ann$Ethnicity))

i = 1
j = 1
k = 1
for (i in 1:length(ICR_clusters)){
  ICR_cluster = ICR_clusters[i]
  for (j in 1:length(IMS_categories)){
    IMS = IMS_categories[j]
    for (k in 1:length(ICR_genes)){
      gene = ICR_genes[k]
      
      plot_df = Clinical_data_ann[, c("Sample_ID","HL.ICR.Cluster", "IMS", "Ethnicity")]
      plot_df$gene = RNASeq.QN.LOG2[gene,][match(plot_df$Sample_ID, colnames(RNASeq.QN.LOG2))]
      plot_df = plot_df[which(plot_df$HL.ICR.Cluster == ICR_cluster & plot_df$IMS == IMS),]
      
      if(nrow(plot_df) > 1){
        my_comparisons = list(c("Asian", "Arab"))
        
        plot = ggplot(plot_df, aes(x = Ethnicity, y = gene)) +
          geom_boxplot(outlier.shape = NA) +
          geom_jitter(width = 0.3, size = 0.5) +
          stat_compare_means(comparisons = my_comparisons, method = "t.test") +
          ylab(gene) +
          theme_bw() +
          ggtitle(paste0(gene, " expression across ethnicities\n",
                         "in ", IMS, " breast cancer in ", ICR_cluster, " cluster"))
        
        png(filename = paste0("./Figures/C11_Boxplot_ICR_gene_across_ethnicities/Boxplot_", ICR_cluster,
                              "_", IMS, "_", gene, "_by_ethnicities.png"), width = 5, height = 5,
            units = "in", res = 600)
        plot(plot)
        dev.off()
      }
    }
  }
}