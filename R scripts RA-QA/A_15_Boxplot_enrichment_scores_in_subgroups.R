
# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-RA-QA")                                                                   # Setwd to location were output files have to be saved.

source("../tools/ipak.function.R")

ipak(c("ggplot2", "ggpubr", "dplyr"))

# Set parameters
Gene_set = "Selected_pathways" #"Bindea_ORIG"  #"Selected_pathways"

# Load data
load(paste0("./Analysis/ssGSEA_scores/", Gene_set, "_ES_scores.Rdata"))
load("./Analysis/Sample_annotations.Rdata")

# Set parameters
var = "[IPA] AMPK Signaling"

# Analysis
Clinical_data_ann = Clinical_data_ann[-which(Clinical_data_ann$Ethnicity == "Other"),]
Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann$HL.ICR.Cluster)),]

Clinical_data_ann$HL.ICR.Cluster = factor(Clinical_data_ann$HL.ICR.Cluster, levels = c("ICR-High", "ICR-Low"))

plot_df = Clinical_data_ann[, c("Sample_ID", "HL.ICR.Cluster", "IMS", "Ethnicity")]
plot_df$var = ES[var, ][match(plot_df$Sample_ID, colnames(ES))]

my_comparisons = list(c("Arab", "Asian"))

plot = ggplot(plot_df, aes(x = Ethnicity, y = var)) +
  geom_boxplot(outlier.shape = NA, color = "grey") +
  geom_jitter(aes(color = Ethnicity), width = 0.15, size = 0.2) +
  scale_color_manual(values = c("#EE82EE", "#7F2D80")) + # Asian:"#7F2D80", Hispanic: "#0080FF"
  facet_grid(HL.ICR.Cluster~IMS, margins = TRUE) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") +
  theme_bw() +
  theme(axis.text.x = element_text(colour = NA),
        strip.background = element_rect(colour="black", fill=NA)) +
  ylab(paste0(var, " enrichment")) +
  ylim(c(0, 0.55))

dir.create(paste0("./Figures/A15_Boxplot_", Gene_set,"_ES_by_ethnicity"), showWarnings = FALSE)
png(filename = paste0("./Figures/A15_Boxplot_", Gene_set, "_ES_by_ethnicity/A15_Boxplot_", var, "_by_ethnicity_Arab&Asian.png"), 
    width = 5, height = 6, units = "in", res = 600)
plot(plot)
dev.off()
