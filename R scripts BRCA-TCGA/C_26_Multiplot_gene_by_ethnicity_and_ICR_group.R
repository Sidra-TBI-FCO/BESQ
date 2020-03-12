
# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA")                                                                   # Setwd to location were output files have to be saved.

source("../tools/ipak.function.R")

ipak(c("ggplot2", "ggpubr", "dplyr", "easyGgplot2"))

# Set parameters
To_plot = c("CAMKK2")

IMS_group = "BasalMyo"

# Load data
load("./Data/expression matrix/BRCA_gene_RNAseq_normalized_Tissue_Type_and_whitelists_filtered.Rdata")
load("./Analysis/Sample_annotations.Rdata")

# Set parameters
Ethnicity_variable = "Assigned_Ethnicity_simplified" #"Ethnicity_reported_simple"
variables = rownames(filtered.norm.RNAseqData)

# Analysis

RNASeq.QN.LOG2 = log(filtered.norm.RNAseqData + 1, 2)
if(Ethnicity_variable == "Ethnicity_reported_simple"){
  Clinical_data_ann = Clinical_data_ann[-which(Clinical_data_ann[, Ethnicity_variable] %in% c("Other", "Asian", "Hispanic")),]
}
if(Ethnicity_variable == "Assigned_Ethnicity_simplified"){
  Clinical_data_ann = Clinical_data_ann[-which(Clinical_data_ann[,Ethnicity_variable] %in% c("Unclear", "Asian")),]
}

Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann$HML.ICR.Cluster)),]
Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann[, Ethnicity_variable])),]
Clinical_data_ann$HML.ICR.Cluster = factor(Clinical_data_ann$HML.ICR.Cluster, levels = c("ICR High", "ICR Medium", "ICR Low"))
Clinical_data_ann = Clinical_data_ann[which(Clinical_data_ann$IMS_Mathews == IMS_group),]

plot_df1 = Clinical_data_ann[, c("bcr_patient_barcode", "HML.ICR.Cluster", "IMS_Mathews")]
plot_df1$Ethnicity = Clinical_data_ann[, Ethnicity_variable]

my_comparisons = list(c("White", "Black"))

i = 1
for(i in 1:length(To_plot)){
  var = To_plot[i]
  plot_df = plot_df1
  plot_df$var = RNASeq.QN.LOG2[var, ][match(plot_df$bcr_patient_barcode, substring(colnames(RNASeq.QN.LOG2), 1, 12))]
  
  min = min(plot_df$var) - abs(min(plot_df$var)*0.1)
  max = max(plot_df$var) + abs(max(plot_df$var)* 0.25)
  distance = max - min
  position = max - (distance*0.1)
  
  px = ggplot(plot_df, aes(x = Ethnicity, y = var)) +
    geom_boxplot(outlier.shape = NA, color = "grey") +
    geom_jitter(aes(color = Ethnicity), width = 0.15, size = 0.2) +
    scale_color_manual(values = c("#FF7E00", "#066C3C")) + # Asian:"#7F2D80", Hispanic: "#0080FF"
    facet_grid(HML.ICR.Cluster~IMS_Mathews, margins = "HML.ICR.Cluster") +
    stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif",
                       label.y = position) +
    theme_bw() +
    theme(axis.text.x = element_text(colour = NA),
          axis.text.y = element_text(colour = "black", size = 7),
          strip.background = element_blank(),
          strip.text = element_blank(),
          legend.position = "none",
          axis.title.y = element_blank()) +
    ggtitle(var) +
    #strip.placement = NULL) +
    ylab(paste0(var, " enrichment")) +
    ylim(c(min, max))
  
  assign(paste0("p", i), px)
  # + geom_text_repel(aes(label = rownames(plot_df)),
  #size = 3.5)
}

plots = paste("p", 1:length(To_plot), sep = "")
list_of_plots = mget(plots)

dir.create(paste0("./Figures/C_26_Multipanel_Boxplot_Facet_individual_genes"), showWarnings = FALSE)
png(paste0("./Figures/C_26_Multipanel_Boxplot_Facet_individual_genes/v3_Multipanel", IMS_group,".png"),
    width = 1, height = 5, units = "in", res = 600)

#width = 8.27, height = 10, units = "in", res = 600)
#png(paste0("./Figures/Multi_panel_boxplots/ggpaired_core_network_analysis_1_", signif.cutoff, "_multiplot_p", 
#number_start,"-p", number_end,".png"), width = 11, height = 11, units = "in", res = 600)
ggplot2.multiplot(plotlist = list_of_plots[1:length(To_plot)], cols=length(To_plot))
dev.off()
