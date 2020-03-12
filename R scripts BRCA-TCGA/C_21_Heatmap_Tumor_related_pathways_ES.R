#######
#
# Heatmap of ES selected pathways
#
#######

# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA")
source("../tools/ipak.function.R")
required.packages <- c("ComplexHeatmap", "impute", "ctc", "amap", "dendextend")
ipak(required.packages)  

# Load data
load("./Analysis/Sample_annotations.Rdata")
load("./Analysis/ssGSEA_scores/Selected_pathways_ES_scores.Rdata")
#ICR_gx_boost = read.csv("./Analysis/Gxboost_trees_analysis_Raghvendra/Results/Feature_Importance_For_ICR.csv", stringsAsFactors = FALSE)
#ICR_gx_boost$Features = paste("[", ICR_gx_boost$Features, sep = "")
#ICR_gx_boost$Features = gsub("\\:", "] ", ICR_gx_boost$Features)
#ICR_gx_boost$Features = gsub("\\_", " ", ICR_gx_boost$Features)
Survival_gx_boost = read.csv("./Analysis/Gxboost_trees_analysis_Raghvendra/Results/Feature_Importance_For_Survival.csv", stringsAsFactors = FALSE)
Survival_gx_boost$Features = paste("[", Survival_gx_boost$Features, sep = "")
Survival_gx_boost$Features = gsub("\\:", "] ", Survival_gx_boost$Features)
Survival_gx_boost$Features = gsub("\\_", " ", Survival_gx_boost$Features)
pathways_1 = c("[TPW] PI3Kgamma Signature", "[TBI] MAPK up genes", "[IPA] UVB Induced MAPK Signaling",
             "[IPA] ErbB2 ErbB3 Signaling", "[IPA] Estrogen Dependent Breast Cancer Signaling", "[IPA] AMPK Signaling",
             "[IPA] PTEN Signaling", "[IPA] ERK MAPK Signaling", "[HM] Estrogen response",
             "[HM] Angiogenesis", "[HM] PI3K Akt mTOR signaling", "[LM] Proliferation",
             "[IPA] Mismatch Repair in Eukaryotes", "[HM] KRAS signaling down", "[HM] UV response up",
             "[HM] Glycolysis", "[HM] E2F targets", "[HM] G2M checkpoint")
pathways_2 = c("[HM] PI3K Akt mTOR signaling", "[IPA] ErbB2 ErbB3 Signaling",
               "[HM] G2M checkpoint", "[IPA] PI3K AKT Signaling", 
               "[HM] KRAS signaling down", "[IPA] ERK MAPK Signaling")

pathways_3 = c("[HM] PI3K Akt mTOR signaling", "[IPA] ErbB2 ErbB3 Signaling",
               "[HM] G2M checkpoint", "[IPA] PI3K AKT Signaling", 
               "[HM] KRAS signaling down", "[IPA] ERK MAPK Signaling",
               "[TBI] MAPK up genes", "[HM] Angiogenesis")

pathways_4 = c("[HM] UV response up", "[HM] Glycolysis", "[HM] KRAS signaling down",
               "[IPA] AMPK Signaling", "[IPA] ErbB2 ErbB3 Signaling")

load("./Analysis/t_test_between_ethnicities/C23_significant_pathways_and_importance_gxboost_tree_survival_BasalMyo_All.Rdata")

# Set parameters
IMS_groups = c("BasalMyo") # c("BasalMyo", "BasalHer2")
Pathways = pathways_4
#Pathways = significant_pathways
Name_pathway_selection = "pathways_4"
#Pathways = Survival_gx_boost$Features[c(1:8)]
#Pathways = c(Pathways, "ICR genes")

# Preparation Heatmap
Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann$IMS_Mathews)),]
Clinical_data_ann = Clinical_data_ann[which(Clinical_data_ann$IMS_Mathews %in% IMS_groups),]
Clinical_data_ann = Clinical_data_ann[which(Clinical_data_ann$Assigned_Ethnicity_simplified %in% c("Black", "White")),]
Clinical_data_ann = Clinical_data_ann[order(Clinical_data_ann$Assigned_Ethnicity_simplified),]

Clinical_data_ann$Ethnicity_ICR = paste(Clinical_data_ann$Assigned_Ethnicity_simplified, Clinical_data_ann$HML.ICR.Cluster)

colnames(ES) = substring(colnames(ES), 1, 12)
ES = ES[,Clinical_data_ann$bcr_patient_barcode]

ES = ES[Pathways,]

ESz = ES
for(j in 1: nrow(ESz))  {
  ESz[j,] = (ES[j,]-mean(ES[j,]))/sd(ES[j,]) # z-score the enrichment matrix
}

Ethnicity_colors = c("White" = "#FFEC42", "Black" = "#066C3C")

#IMS_Mathews_colors = c("BasalHer2" = "#F4A101", "BasalLumHer2" = "#92D050",
#                       "BasalMyo" = "#3D64E4", "MyoLumB" = "#C720F1",
#                       "Lum" = "#158400","LumBasal" = "#D36733",
#                       "MyoLumA" = "#6F1287", "MyoLumHer2" = "#F377C6")

Ethnicity_Black_ICR_colors = c("Black ICR High" = "red", "Black ICR Low" = "blue", "Black ICR Medium" = "green",
                               "White ICR High" = "grey", "White ICR Low" = "grey", "White ICR Medium" = "grey")

Ethnicity_White_ICR_colors = c("Black ICR High" = "grey", "Black ICR Low" = "grey", "Black ICR Medium" = "grey",
                               "White ICR High" = "red", "White ICR Low" = "blue", "White ICR Medium" = "green")

#col_fun = circlize::colorRamp2(c(min(Expression.matrix.z), max(Expression.matrix.z)), c("blue", "red"))
col_fun = circlize::colorRamp2(c(-4,-0.1,0,0.1, 3), c("blue", "#EAE7EF","white","#FAD9D0", "red"))

mat =  ESz
#ha_column = HeatmapAnnotation(df = data.frame(IMS_Mathews = Clinical_data_ann$IMS_Mathews,
#                                              Ethnicity = Clinical_data_ann$Assigned_Ethnicity_simplified),
#                              show_annotation_name = TRUE,
#                              col = list(IMS_Mathews = IMS_Mathews_colors,
#                                         Ethnicity = Ethnicity_colors))

ha_column = HeatmapAnnotation(df = data.frame(Ethnicity = Clinical_data_ann$Assigned_Ethnicity_simplified,
                                              Ethnicity_ICR_Black = Clinical_data_ann$Ethnicity_ICR,
                                              Ethnicity_ICR_White = Clinical_data_ann$Ethnicity_ICR),
                              show_annotation_name = TRUE,
                              col = list(Ethnicity = Ethnicity_colors,
                                         Ethnicity_ICR_Black = Ethnicity_Black_ICR_colors,
                                         Ethnicity_ICR_White = Ethnicity_White_ICR_colors))

## Cluster tree
col_dend = as.dendrogram(hclust(dist(t(ESz))), method = "complete")
#col_dend = color_branches(col_dend, k = 2, col = c( "purple", "orange"))
x = as.data.frame(t(ESz))
x$cluster = cutree(col_dend, k=2)[match(rownames(x), names(cutree(col_dend, k=2)))]
x$cluster[which(x$cluster == 1)] = "1"
x$cluster[which(x$cluster == 2)] = "2"

x$Ethnicity = Clinical_data_ann$Assigned_Ethnicity_simplified[match(rownames(x),
                                                                    Clinical_data_ann$bcr_patient_barcode)]

dir.create("./Analysis/Benefit_Cluster", showWarnings = FALSE)
save(x, file = paste0("./Analysis/Benefit_Cluster/", Name_pathway_selection , "_", paste0(IMS_groups, collapse = "_"), ".Rdata"))

# Heatmap
#png(paste0("./Figures/Heatmaps/",  Gene.set.list, "_", Gene.set,".png"), res = 600, width = 10, height = 10, units = "in")

dir.create("./Figures/Heatmaps/C_21_Heatmap_Selected_Pathways_ES/", showWarnings = FALSE)
png(paste0("./Figures/Heatmaps/C_21_Heatmap_Selected_Pathways_ES/C_21_", Name_pathway_selection,"_Heatmap_Selected_Pathways_ES", 
           paste0(IMS_groups, collapse = "_"), "v1.png"), res = 600, width = 7, height = 4, units = "in")
plot = Heatmap(mat,
               name = "expression values", 
               show_heatmap_legend = FALSE,
               cluster_rows = TRUE,
               cluster_columns = col_dend,
               #cluster_columns = FALSE,
               row_title_gp = gpar(fontsize = 0.1),
               #column_names_gp = gpar(fontsize = 10),
               row_names_gp = gpar(fontsize = 10),
               top_annotation = ha_column,
               #column_title = paste0("Expression PAM50 genes"),
               show_column_names = FALSE,
               show_row_names = TRUE
) 

draw(plot, show_annotation_legend = TRUE, gap = unit(1, "mm"))

dev.off()




