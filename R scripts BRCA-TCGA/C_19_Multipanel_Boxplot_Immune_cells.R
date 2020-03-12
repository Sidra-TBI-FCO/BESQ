
# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA")                                                                   # Setwd to location were output files have to be saved.

source("../tools/ipak.function.R")

ipak(c("ggplot2", "ggpubr", "dplyr", "easyGgplot2"))

# Set parameters
Gene_set = "Bindea_ORIG" #"Bindea_ORIG"  #"Selected_pathways"
Combine_ICR_Medium_Low = "Combine_ICR_Medium_Low"
load("./Analysis/t_test_between_ethnicities/C23_significant_pathways_BasalMyo_All.Rdata")

To_plot = c("DC", "TReg", "Th2 cells", "B cells")
  # c("[HM] UV response up", "[HM] Glycolysis", "[HM] KRAS signaling down",
  # "[IPA] AMPK Signaling", "[IPA] ErbB2 ErbB3 Signaling")
  #c("[HM] G2M checkpoint", "[HM] KRAS signaling down",
   #         "[TBI] MAPK up genes", "[HM] Angiogenesis",
    #        "[HM] PI3K Akt mTOR signaling", "[IPA] ErbB2 ErbB3 Signaling",
     #       "[IPA] ERK MAPK Signaling", "[IPA] PI3K AKT Signaling") #c("MDSC", "Activated CD8 T cell", "Activated CD4 T cell" )
  
  
  # c("[HM] UV response up", "[HM] Glycolysis", "[HM] KRAS signaling down",
            # "[IPA] AMPK Signaling", "[IPA] ErbB2 ErbB3 Signaling")
# c("[HM] PI3K Akt mTOR signaling",
#  "[IPA] ErbB Signaling",
# "[IPA] ErbB2 ErbB3 Signaling", "[IPA] AMPK Signaling",
# "[HM] Estrogen response", "[TBI] MAPK up genes",
# "[HM] Angiogenesis", "[IPA] Estrogen Dependent Breast Cancer Signaling") 
#c("B cells", "DC", "NK CD56dim cells", "TReg", "Th2 cells", "NK cells")
IMS_group = "BasalMyo"

# Load data
load(paste0("./Analysis/ssGSEA_scores/", Gene_set, "_ES_scores.Rdata"))
load("./Analysis/Sample_annotations.Rdata")
load(paste0("./Analysis/t_test_between_ethnicities/C14_Diff_", Gene_set, "_between_ICR_IMS_Ethnicity.Rdata"))

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

if(Combine_ICR_Medium_Low == "Combine_ICR_Medium_Low"){
  levels(Clinical_data_ann$HML.ICR.Cluster) = c("ICR High", "ICR Medium-Low", "ICR Medium-Low")
}
Clinical_data_ann = Clinical_data_ann[which(Clinical_data_ann$IMS_Mathews == IMS_group),]

plot_df1 = Clinical_data_ann[, c("bcr_patient_barcode", "HML.ICR.Cluster", "IMS_Mathews")]
plot_df1$Ethnicity = Clinical_data_ann[, Ethnicity_variable]

my_comparisons = list(c("White", "Black"))

results_df = data.frame(Celltype = To_plot, p_value = NA)
results_df1 = data.frame(Celltype = To_plot, p_value = NA)
results_df3 = data.frame(Celltype = To_plot, p_value = NA)

i = 1
for(i in 1:length(To_plot)){
  var = To_plot[i]
  plot_df = plot_df1
  plot_df$var = ES[var, ][match(plot_df$bcr_patient_barcode, substring(colnames(ES), 1, 12))]
  
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
                       label.y = position, hide.ns = FALSE) +
    theme_bw() +
    theme(axis.text.x = element_text(colour = NA),
          axis.text.y = element_text(colour = "black", size = 7),
          strip.background = element_blank(),
          strip.text = element_blank(),
          legend.position = "none",
          axis.title.y = element_blank(),
          aspect.ratio =1.6/1) +
    ggtitle(var) +
          #strip.placement = NULL) +
    ylab(paste0(var, " enrichment")) +
    ylim(c(min, max)) +
    xlab("")
  
  assign(paste0("p", i), px)
  # + geom_text_repel(aes(label = rownames(plot_df)),
  #size = 3.5)
  plot_df2 = plot_df[which(plot_df$HML.ICR.Cluster == "ICR Medium-Low"),]
  plot_df3 = plot_df[which(plot_df$HML.ICR.Cluster == "ICR High"),]
  
  t_test = t.test(plot_df2$var[which(plot_df2$Ethnicity == "White")], plot_df2$var[which(plot_df2$Ethnicity == "Black")])
  results_df$p_value[which(results_df$Celltype == var)] = t_test$p.value
  
  t_test1 = t.test(plot_df$var[which(plot_df$Ethnicity == "White")], plot_df$var[which(plot_df$Ethnicity == "Black")])
  results_df1$p_value[which(results_df1$Celltype == var)] = t_test1$p.value
  
  t_test3 = t.test(plot_df3$var[which(plot_df3$Ethnicity == "White")], plot_df3$var[which(plot_df3$Ethnicity == "Black")])
  results_df3$p_value[which(results_df3$Celltype == var)] = t_test3$p.value
}

results_df$FDR = p.adjust(p = results_df$p_value, method = "fdr",n = 24*3)
results_df1$FDR = p.adjust(p = results_df1$p_value, method = "fdr",n = 24*3)
results_df3$FDR = p.adjust(p = results_df1$p_value, method = "fdr",n = 24*3)

plots = paste("p", 1:length(To_plot), sep = "")
list_of_plots = mget(plots)

dir.create(paste0("./Figures/C_19_Multipanel_Boxplot_Facet_", Gene_set), showWarnings = FALSE)
png(paste0("./Figures/C_19_Multipanel_Boxplot_Facet_", Gene_set, "/2_oct_2019_Figure3_Final_Multipanel", 
           IMS_group, "_", Combine_ICR_Medium_Low,".png"),
    width = 5, height = 4, units = "in", res = 600)
    
    #width = 8.27, height = 10, units = "in", res = 600)
#png(paste0("./Figures/Multi_panel_boxplots/ggpaired_core_network_analysis_1_", signif.cutoff, "_multiplot_p", 
#number_start,"-p", number_end,".png"), width = 11, height = 11, units = "in", res = 600)
ggplot2.multiplot(plotlist = list_of_plots[1:length(To_plot)], cols=length(To_plot))
dev.off()
