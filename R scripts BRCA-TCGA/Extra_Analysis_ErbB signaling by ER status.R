
# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA")                                                                   # Setwd to location were output files have to be saved.

source("../tools/ipak.function.R")

ipak(c("ggplot2", "ggpubr", "dplyr"))

# Set parameters
Gene_set = "Selected.pathways" #"Bindea_ORIG"  #"Selected_pathways"
Stages = ""

# Load data
load(paste0("./Analysis/ssGSEA_scores/v2_", Gene_set, "_ES_scores.Rdata"))
load("./Analysis/Sample_annotations.Rdata")
load("./Analysis/t_test_between_ethnicities/C23_significant_pathways_BasalMyo_All.Rdata")
patient_data = read.csv("./Data/Clinical Data/patient.csv", stringsAsFactors = FALSE)
patient_data = patient_data[-c(1, 2),]
#load(paste0("./Analysis/t_test_between_ethnicities/C14_Diff_", Gene_set, "_between_ICR_IMS_Ethnicity.Rdata"))

# Set parameters
Ethnicity_variable = "Assigned_Ethnicity_simplified" #"Ethnicity_reported_simple"
#variables = significant_pathways
variables = rownames(ES)

# Analysis
if(Ethnicity_variable == "Ethnicity_reported_simple"){
  Clinical_data_ann = Clinical_data_ann[-which(Clinical_data_ann[, Ethnicity_variable] %in% c("Other", "Asian", "Hispanic")),]
}
if(Ethnicity_variable == "Assigned_Ethnicity_simplified"){
  Clinical_data_ann = Clinical_data_ann[-which(Clinical_data_ann[,Ethnicity_variable] %in% c("Unclear", "Asian")),]
}

#Clinical_data_ann = Clinical_data_ann[which(Clinical_data_ann$IMS_Mathews == "BasalMyo"),]
table(Clinical_data_ann$IMS_Mathews, exclude = NULL)
Clinical_data_ann$IMS_Mathews = as.character(Clinical_data_ann$IMS_Mathews)
Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann$IMS_Mathews)),]
#Clinical_data_ann$IMS_Mathews[-which(Clinical_data_ann$IMS_Mathews == "BasalMyo")] = "Non-BasalMyo"
table(Clinical_data_ann$IMS_Mathews, exclude = NULL)

if(Stages == c("StageI_II")){
  Clinical_data_ann = Clinical_data_ann[which(Clinical_data_ann$ajcc_pathologic_tumor_stage %in% c("Stage I", "Stage IA",
                                                                                                   "Stage II", "Stage IIA", "Stage IIB")),]
}

#Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann$HML.ICR.Cluster)),]
Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann[, Ethnicity_variable])),]

Clinical_data_ann$HML.ICR.Cluster = factor(Clinical_data_ann$HML.ICR.Cluster, levels = c("ICR High", "ICR Medium", "ICR Low"))

Clinical_data_ann$Estrogen_receptor = patient_data$er_status_by_ihc[match(Clinical_data_ann$bcr_patient_barcode,
                                                                           patient_data$bcr_patient_barcode)]

table(Clinical_data_ann$Estrogen_receptor, exclude = NULL)
Clinical_data_ann$Estrogen_receptor[which(Clinical_data_ann$Estrogen_receptor %in% c("[Not Evaluated]", "Indeterminate"))] = NA
#Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann$Estrogen_receptor)),]
Clinical_data_ann$Category = NA
Clinical_data_ann$Category[which(Clinical_data_ann$Estrogen_receptor == "Positive")] = "ER positive"
Clinical_data_ann$Category[which(Clinical_data_ann$IMS_Mathews == "BasalMyo")] = "BasalMyo"

Clinical_data_ann = Clinical_data_ann[which(Clinical_data_ann$IMS_Mathews == "BasalMyo"),]

i = 56 # 2, 5, 7, 14
for (i in 1:length(variables)){
  var = variables[i]
  plot_df = Clinical_data_ann[, c("bcr_patient_barcode", "IMS_Mathews", "Estrogen_receptor")]
  plot_df$Ethnicity = Clinical_data_ann[, Ethnicity_variable]
  plot_df$var = ES[var, ][match(plot_df$bcr_patient_barcode, substring(colnames(ES), 1, 12))]
  
  my_comparisons = list(c("White", "Black"))
  #my_comparisons = list(c("Positive", "Negative"))
  table(plot_df$Estrogen_receptor) 
  
  min = min(plot_df$var) - abs(min(plot_df$var)*0.1)
  max = max(plot_df$var) + abs(max(plot_df$var)* 0.08)
  max = max(plot_df$var) + abs(max(plot_df$var)* 0.09)
  plot = ggplot(plot_df, aes(x = Ethnicity, y = var)) +
    geom_boxplot(outlier.shape = NA, color = "grey") +
    geom_jitter(aes(color = Ethnicity), width = 0.15, size = 0.2) +
    scale_color_manual(values = c("#FF7E00", "#066C3C")) + # Asian:"#7F2D80", Hispanic: "#0080FF"
    stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") +
    theme_bw() +
    theme(axis.text.x = element_text(colour = "black", size =10, angle = 45, vjust = 1, hjust = 1),
          axis.text.y = element_text(colour = "black", size = 15),
          axis.title.y = element_text(colour = "black", size = 15),
          axis.title.x = element_text(colour = "black", size = 15),
          strip.background = element_rect(colour="black", fill=NA),
          legend.position = "none",
          aspect.ratio = 3/1.5) +
    ylab(paste0(var)) +
    xlab("ER status") +
    ylim(c(min, max)) #+
    #facet_grid(.~IMS_Mathews)
  
  dir.create(paste0("./Figures/Reviewer_comment_", Gene_set,"_ES_by_ethnicity"), showWarnings = FALSE)
  png(filename = paste0("./Figures/Reviewer_comment_", Gene_set, "_ES_by_ethnicity/v5_", gsub("./", "_", var), "_by_", Ethnicity_variable,"_and_ER_positive_Black&White.png"), 
      width = 2.5, height = 3.2, units = "in", res = 600)
  plot(plot)
  dev.off()
}

i=5
for (i in 1:length(variables)){
  var = variables[i]
  plot_df = Clinical_data_ann[, c("bcr_patient_barcode", "Category")]
  plot_df$Ethnicity = Clinical_data_ann[, Ethnicity_variable]
  plot_df = plot_df[-which(is.na(plot_df$Category)),]
  
  plot_df$var = ES[var, ][match(plot_df$bcr_patient_barcode, substring(colnames(ES), 1, 12))]
  
  my_comparisons = list(c("White", "Black"))
  table(plot_df$Category) 
  
  min = min(plot_df$var) - abs(min(plot_df$var)*0.1)
  max = max(plot_df$var) + abs(max(plot_df$var)* 0.08)
  plot = ggplot(plot_df, aes(x = Category, y = var)) +
    geom_boxplot(outlier.shape = NA, color = "grey") +
    geom_jitter(aes(color = Ethnicity), width = 0.15, size = 0.2) +
    scale_color_manual(values = c("#FF7E00", "#066C3C")) + # Asian:"#7F2D80", Hispanic: "#0080FF"
    stat_compare_means(comparisons = list(c("BasalMyo", "ER positive")), method = "t.test", label = "p.signif") +
    theme_bw() +
    theme(axis.text.x = element_text(colour = "black", size =12, angle = 45, vjust = 1, hjust = 1),
          axis.text.y = element_text(colour = "black", size = 12),
          axis.title.y = element_text(colour = "black", size = 12),
          axis.title.x = element_text(colour = "black", size = 12),
          strip.background = element_rect(colour="black", fill=NA),
          legend.position = "none",
          aspect.ratio = 3/1.5) +
    ylab(paste0(var)) +
    xlab("") +
    ylim(c(min, max)) 
  
  dir.create(paste0("./Figures/Reviewer_comment_", Gene_set,"_ES_by_ethnicity"), showWarnings = FALSE)
  png(filename = paste0("./Figures/Reviewer_comment_", Gene_set, "_ES_by_ethnicity/v2_", gsub("./", "_", var), "_by_", Ethnicity_variable,"_and_ER_positive_Black&White.png"), 
      width = 2.5, height = 3.2, units = "in", res = 600)
  plot(plot)
  dev.off()
}

