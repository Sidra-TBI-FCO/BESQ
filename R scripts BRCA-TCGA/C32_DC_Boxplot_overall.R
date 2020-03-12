
# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA")                                                                   # Setwd to location were output files have to be saved.

source("../tools/ipak.function.R")

ipak(c("ggplot2", "ggpubr", "dplyr"))

# Set parameters
Gene_set = "Selected_pathways" #"Bindea_ORIG"  #"Selected_pathways"

# Load data
load(paste0("./Analysis/ssGSEA_scores/", Gene_set, "_ES_scores.Rdata"))
load("./Analysis/Sample_annotations.Rdata")

# Set parameters
color_basal_myo = ""
var = "[IPA] AMPK Signaling"
#var = "TReg"  # "TReg" or "Th2 cells"
only_nonBasalMyo = "only_nonBasalMyo"

# Analysis
Clinical_data_ann = Clinical_data_ann[which(Clinical_data_ann$Assigned_Ethnicity_simplified %in% c("White", "Black", "Asian")),]
#Clinical_data_ann$Ethnicity[which(Clinical_data_ann$Ethnicity %in% c("Asian", "Other"))] = "Other"
#Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann$HL.ICR.Cluster)),]

#Clinical_data_ann$HL.ICR.Cluster = factor(Clinical_data_ann$HL.ICR.Cluster, levels = c("ICR-High", "ICR-Low"))
if(only_nonBasalMyo == "only_nonBasalMyo"){
  Clinical_data_ann = Clinical_data_ann[-which(Clinical_data_ann$IMS_Mathews == "BasalMyo"),]
}

plot_df = Clinical_data_ann[,c("bcr_patient_barcode", "Assigned_Ethnicity_simplified")]
colnames(plot_df) = c("Sample_ID", "Ethnicity")
plot_df$var = ES[var, ][match(plot_df$Sample_ID, substring(colnames(ES), 1, 12))]
plot_df$IMS_Mathews = Clinical_data_ann$IMS_Mathews[match(plot_df$Sample_ID, Clinical_data_ann$bcr_patient_barcode)]

if(color_basal_myo == "color_basal_myo"){
  plot_df$Othnicity = plot_df$Ethnicity
  plot_df$Othnicity = as.character(plot_df$Othnicity)
  plot_df$IMS_Mathews = as.character(plot_df$IMS_Mathews)
  plot_df$Othnicity[which(plot_df$IMS_Mathews == "BasalMyo")] = "B"
  plot_df$Othnicity = factor(plot_df$Othnicity, levels = c("White", "Black", "Asian", "B"))
  colors = c("#FF7E00", "#066C3C", "#7F2D80", "green")
}
if(color_basal_myo == ""){
  colors = c("#FF7E00", "#066C3C", "#7F2D80")
}

my_comparisons = list(c("Black", "White"), c("Asian", "White"), c("Black", "Asian"))

plot = ggplot(plot_df, aes(x = Ethnicity, y = var)) +
  geom_boxplot(outlier.shape = NA, color = "grey") +
  geom_jitter(aes(color = Ethnicity), width = 0.15, size = 0.1) +
  scale_color_manual(values = colors) + # Asian:"#7F2D80", Hispanic: "#0080FF"
  #stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") +
  theme_bw() +
  #facet_grid(.~IMS_Mathews) +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        strip.background = element_rect(colour="black", fill=NA)) +
  ylab(paste0(var, " enrichment")) + 
  ylim(0.14, 0.26)

dir.create(paste0("./Figures/C32_Boxplot_", Gene_set,"_ES_by_ethnicity"), showWarnings = FALSE)
png(filename = paste0("./Figures/C32_Boxplot_", Gene_set, "_ES_by_ethnicity/", only_nonBasalMyo, "_", color_basal_myo, "_Oct_C32_v3_Boxplot_", var, "_by_ethnicity_Arab&Asian.png"), 
    width = 2.5, height = 1.9, units = "in", res = 600)
plot(plot)
dev.off()

# facet bigger plot
png(filename = paste0("./Figures/C32_Boxplot_", Gene_set, "_ES_by_ethnicity/", only_nonBasalMyo, "_", color_basal_myo, "_C32_v4_facet_Boxplot_", var, "_by_ethnicity_Arab&Asian.png"), 
    width = 8, height = 5, units = "in", res = 600)
plot(plot)
dev.off()
