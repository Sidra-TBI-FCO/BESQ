
# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA")                                                                   # Setwd to location were output files have to be saved.

source("../tools/ipak.function.R")

ipak(c("ggplot2", "ggpubr", "dplyr"))

# Set parameters
Gene_set = "Selected_pathways" #"Bindea_ORIG"  #"Selected_pathways"
Stages = "StageI_II"

# Load data
load(paste0("./Analysis/ssGSEA_scores/", Gene_set, "_ES_scores.Rdata"))
load("./Analysis/Sample_annotations.Rdata")
load("./Analysis/t_test_between_ethnicities/C23_significant_pathways_BasalMyo_All.Rdata")
#load(paste0("./Analysis/t_test_between_ethnicities/C14_Diff_", Gene_set, "_between_ICR_IMS_Ethnicity.Rdata"))

# Set parameters
Ethnicity_variable = "Assigned_Ethnicity_simplified" #"Ethnicity_reported_simple"
variables = significant_pathways
#variables = rownames(ES)

# Analysis
if(Ethnicity_variable == "Ethnicity_reported_simple"){
  Clinical_data_ann = Clinical_data_ann[-which(Clinical_data_ann[, Ethnicity_variable] %in% c("Other", "Asian", "Hispanic")),]
}
if(Ethnicity_variable == "Assigned_Ethnicity_simplified"){
  Clinical_data_ann = Clinical_data_ann[-which(Clinical_data_ann[,Ethnicity_variable] %in% c("Unclear", "Asian")),]
}

Clinical_data_ann = Clinical_data_ann[which(Clinical_data_ann$IMS_Mathews == "BasalMyo"),]

if(Stages == c("StageI_II")){
  Clinical_data_ann = Clinical_data_ann[which(Clinical_data_ann$ajcc_pathologic_tumor_stage %in% c("Stage I", "Stage IA",
                                                                                             "Stage II", "Stage IIA", "Stage IIB")),]
}

#Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann$HML.ICR.Cluster)),]
Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann[, Ethnicity_variable])),]

Clinical_data_ann$HML.ICR.Cluster = factor(Clinical_data_ann$HML.ICR.Cluster, levels = c("ICR High", "ICR Medium", "ICR Low"))

i = 11
for (i in 1:length(variables)){
  var = variables[i]
  plot_df = Clinical_data_ann[, c("bcr_patient_barcode", "IMS_Mathews")]
  plot_df$Ethnicity = Clinical_data_ann[, Ethnicity_variable]
  plot_df$var = ES[var, ][match(plot_df$bcr_patient_barcode, substring(colnames(ES), 1, 12))]
  
  my_comparisons = list(c("White", "Black"))
  
  min = min(plot_df$var) - abs(min(plot_df$var)*0.1)
  max = max(plot_df$var) + abs(max(plot_df$var)* 0.08)
  plot = ggplot(plot_df, aes(x = Ethnicity, y = var)) +
    geom_boxplot(outlier.shape = NA, color = "grey") +
    geom_jitter(aes(color = Ethnicity), width = 0.15, size = 0.2) +
    scale_color_manual(values = c("#FF7E00", "#066C3C")) + # Asian:"#7F2D80", Hispanic: "#0080FF"
    stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") +
    theme_bw() +
    theme(axis.text.x = element_text(colour = NA),
          axis.text.y = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black", size = 5),
          strip.background = element_rect(colour="black", fill=NA),
          legend.position = "none",
          aspect.ratio = 3/1.5) +
    ylab("") +
    xlab("") +
    ylim(c(min, max))
  
  dir.create(paste0("./Figures/C15_Boxplot_", Gene_set,"_ES_by_ethnicity"), showWarnings = FALSE)
  png(filename = paste0("./Figures/C15_Boxplot_", Gene_set, "_ES_by_ethnicity/", Stages, "_C15.1_Boxplot_", gsub("./", "_", var), "_by_", Ethnicity_variable,"_Black&White.png"), 
      width = 1.5, height = 2, units = "in", res = 600)
  plot(plot)
  dev.off()
}
