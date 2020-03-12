
# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA")                                                                   # Setwd to location were output files have to be saved.

source("../tools/ipak.function.R")

ipak(c("ggplot2", "ggpubr", "dplyr"))

# Set parameters
Gene_set = "Bindea_ORIG" #"Bindea_ORIG"  #"Selected_pathways"

# Load data
load(paste0("./Analysis/ssGSEA_scores/", Gene_set, "_ES_scores.Rdata"))
load("./Analysis/Sample_annotations.Rdata")
#load(paste0("./Analysis/t_test_between_ethnicities/C14_Diff_", Gene_set, "_between_ICR_IMS_Ethnicity.Rdata"))

# Set parameters
Ethnicity_variable = "Assigned_Ethnicity_simplified" #"Ethnicity_reported_simple"
variables = rownames(ES)

# Analysis
if(Ethnicity_variable == "Ethnicity_reported_simple"){
  Clinical_data_ann = Clinical_data_ann[-which(Clinical_data_ann[, Ethnicity_variable] %in% c("Other", "Asian", "Hispanic")),]
}
if(Ethnicity_variable == "Assigned_Ethnicity_simplified"){
  Clinical_data_ann = Clinical_data_ann[-which(Clinical_data_ann[,Ethnicity_variable] %in% c("Unclear", "Asian")),]
}

Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann$HML.ICR.Cluster)),]
Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann[, Ethnicity_variable])),]

Clinical_data_ann$HML.ICR.Cluster = factor(Clinical_data_ann$HML.ICR.Cluster, levels = c("ICR High", "ICR Medium", "ICR Low"))

results_df = data.frame(Celltypes = variables, p_value = NA)
results_df_BasalMyo = data.frame(Celltypes = variables, p_value = NA)
results_df_BasalMyo_ML = data.frame(Celltypes = variables, p_value = NA)

i = 15
for (i in 1:length(variables)){
  var = variables[i]
  plot_df = Clinical_data_ann[, c("bcr_patient_barcode", "HML.ICR.Cluster", "IMS_Mathews")]
  plot_df$Ethnicity = Clinical_data_ann[, Ethnicity_variable]
  plot_df$var = ES[var, ][match(plot_df$bcr_patient_barcode, substring(colnames(ES), 1, 12))]
  
  t_test = t.test(plot_df$var[which(plot_df$Ethnicity == "White")], plot_df$var[which(plot_df$Ethnicity == "Black")])
  results_df$p_value[which(results_df$Celltypes == var)] = t_test$p.value
  
  basalmyo_df = plot_df[which(plot_df$IMS_Mathews == "BasalMyo"),]
  t_test2 = t.test(basalmyo_df$var[which(basalmyo_df$Ethnicity == "White")], basalmyo_df$var[which(basalmyo_df$Ethnicity == "Black")])
  results_df_BasalMyo$p_value[which(results_df_BasalMyo$Celltypes == var)] = t_test2$p.value
  
  my_comparisons = list(c("White", "Black"))
  
  min = min(plot_df$var) - abs(min(plot_df$var)*0.1)
  max = max(plot_df$var) + abs(max(plot_df$var)* 0.25)
  plot = ggplot(plot_df, aes(x = Ethnicity, y = var)) +
    geom_boxplot(outlier.shape = NA, color = "grey") +
    geom_jitter(aes(color = Ethnicity), width = 0.15, size = 0.2) +
    scale_color_manual(values = c("#FF7E00", "#066C3C")) + # Asian:"#7F2D80", Hispanic: "#0080FF"
    facet_grid(HML.ICR.Cluster~IMS_Mathews, margins = TRUE) +
    stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") +
    theme_bw() +
    theme(axis.text.x = element_text(colour = NA),
          strip.background = element_rect(colour="black", fill=NA)) +
    ylab(paste0(var, " enrichment")) +
    ylim(c(min, max))
  
  dir.create(paste0("./Figures/C15_Boxplot_", Gene_set,"_ES_by_ethnicity"), showWarnings = FALSE)
  #png(filename = paste0("./Figures/C15_Boxplot_", Gene_set, "_ES_by_ethnicity/C15_v3_Boxplot_", gsub("./", "_", var), "_by_", Ethnicity_variable,"_Black&White.png"), 
      #width = 6.5, height = 6.5, units = "in", res = 600)
  #plot(plot)
  #dev.off()
  
}

