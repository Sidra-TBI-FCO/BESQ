

# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA")                                                                   # Setwd to location were output files have to be saved.

source("../tools/ipak.function.R")

ipak(c("ggplot2", "ggpubr"))

# Set parameters
Ethnicity_variable = "Assigned_Ethnicity_simplified" #"Ethnicity_reported_simple"
ICR_variable = "ICR_ES"

# Load data
load("./Analysis/Sample_annotations.Rdata")
#load("./Analysis/ICR data/TCGA_BRCA_table_cluster_assignment.RData")

#Clinical_data_ann$ICRscore = clustering$ICRscore[match(Clinical_data_ann$bcr_patient_barcode,
                                                       #substring(rownames(clustering), 1, 12))]
#Clinical_data_ann$HML.ICR.Cluster = clustering$HML.ICR.Cluster[match(Clinical_data_ann$bcr_patient_barcode,
                                                             # substring(rownames(clustering), 1, 12))]
#save(Clinical_data_ann, file = "./Analysis/Sample_annotations.Rdata")


Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann[, Ethnicity_variable])),]
if(Ethnicity_variable == "Ethnicity_reported_simple"){
  Clinical_data_ann = Clinical_data_ann[-which(Clinical_data_ann[, Ethnicity_variable] == "Other"),]
}
if(Ethnicity_variable == "Assigned_Ethnicity_simplified"){
  Clinical_data_ann = Clinical_data_ann[-which(Clinical_data_ann[, Ethnicity_variable] == "Unclear"),]
}

Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann$ICRscore)),]

plot_df = data.frame(Patient_id = Clinical_data_ann$bcr_patient_barcode,
                     Ethnicity = Clinical_data_ann[, Ethnicity_variable],
                     IMS = Clinical_data_ann$IMS_Mathews,
                     ICR = Clinical_data_ann[,ICR_variable])

my_comparisons = list(c("White", "Black"), c("White", "Asian"),
                      c("Black", "Asian"))

plot = ggplot(plot_df, aes(y = ICR, x = Ethnicity)) +
  geom_boxplot(outlier.shape = NA, color = "grey") +
  geom_jitter(aes(color = Ethnicity), width = 0.1, size = 0.5) +
  scale_color_manual(values = c("#FF7E00", "#066C3C",  "#7F2D80", "#0080FF")) +
  #stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") +
  facet_grid(.~IMS) +
  theme_bw() +
  ylim(c(-0.3, 0.35)) +
  theme(axis.text.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = NA),
        strip.background = element_rect(colour="black", fill=NA))


dir.create("./Figures/C9_Boxplot_ICRscore_by_ethnicity_facet_for_ethnicity", showWarnings = FALSE)
png(paste0("./Figures/C9_Boxplot_ICRscore_by_ethnicity_facet_for_ethnicity/C9_Boxplot_stats_v2_", ICR_variable,"_per_", 
           Ethnicity_variable, "_by_IMS.png"), width = 6, height = 3, res = 600, units = "in")
print(plot)
dev.off()
