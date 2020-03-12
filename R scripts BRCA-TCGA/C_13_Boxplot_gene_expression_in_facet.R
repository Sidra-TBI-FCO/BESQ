
# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA")                                                                   # Setwd to location were output files have to be saved.

source("../tools/ipak.function.R")

ipak(c("ggplot2", "ggpubr", "dplyr"))

# Set parameters
Ethnicity_variable = "Assigned_Ethnicity_simplified" #"Ethnicity_reported_simple"
gene = "PDCD1"

# Load data
load("./Data/expression matrix/BRCA_gene_RNAseq_normalized_Tissue_Type_and_whitelists_filtered.Rdata")
load("./Analysis/Sample_annotations.Rdata")
load("./Analysis/t_test_between_ethnicities/C11_Diff_ICR_expression_between_ICR_IMS_Ethnicity.Rdata")

# Analysis
RNAseq_log2 = log(filtered.norm.RNAseqData +1, 2)

Clinical_data_ann = Clinical_data_ann[-which(Clinical_data_ann[, Ethnicity_variable] %in% c("Unclear", "Asian")),]
Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann$HML.ICR.Cluster)),]
Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann[, Ethnicity_variable])),]

Clinical_data_ann$HML.ICR.Cluster = factor(Clinical_data_ann$HML.ICR.Cluster, levels = c("ICR High", "ICR Medium", "ICR Low"))

plot_df = Clinical_data_ann[, c("bcr_patient_barcode", "HML.ICR.Cluster", "IMS_Mathews")]
plot_df$Ethnicity = Clinical_data_ann[, Ethnicity_variable]
plot_df$gene = RNAseq_log2[gene, ][match(plot_df$bcr_patient_barcode, substring(colnames(RNAseq_log2), 1, 12))]

my_comparisons = list(c("White", "Black"))

plot = ggplot(plot_df, aes(x = Ethnicity, y = gene)) +
  geom_boxplot(outlier.shape = NA, color = "grey") +
  geom_jitter(aes(color = Ethnicity), width = 0.15, size = 0.2) +
  scale_color_manual(values = c("#FF7E00", "#066C3C")) + # Asian:"#7F2D80", Hispanic: "#0080FF"
  facet_grid(HML.ICR.Cluster~IMS_Mathews, margins = TRUE) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") +
  theme_bw() +
  theme(axis.text.x = element_text(colour = NA),
        strip.background = element_rect(colour="black", fill=NA)) +
  ylab(paste0(gene, " expression")) +
  ylim(c(0, 17))

dir.create("./Figures/C13_Boxplot_ICR_gene_by_ethnicity", showWarnings = FALSE)
png(filename = paste0("./Figures/C13_Boxplot_ICR_gene_by_ethnicity/C13_", Ethnicity_variable, "_Boxplot_", gene, "_by_ethnicity_Black&White.png"), 
    width = 6.5, height = 6.5, units = "in", res = 600)
plot(plot)
dev.off()
