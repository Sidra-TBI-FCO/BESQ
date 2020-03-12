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

# Set parameters
IMS_groups = c("Basal") # c("BasalMyo", "BasalHer2")
# Pathways = c("[TPW] PI3Kgamma Signature", "[TBI] MAPK up genes", "[IPA] UVB Induced MAPK Signaling",
#"[IPA] ErbB2 ErbB3 Signaling", "[IPA] Estrogen Dependent Breast Cancer Signaling", "[IPA] AMPK Signaling",
#"[IPA] PTEN Signaling", "[IPA] ERK MAPK Signaling", "[HM] Estrogen response",
#"[HM] Angiogenesis", "[HM] PI3K Akt mTOR signaling", "[LM] Proliferation",
#"[IPA] Mismatch Repair in Eukaryotes", "[HM] KRAS signaling down", "[HM] UV response up",
#"[HM] Glycolysis", "[HM] E2F targets", "[HM] G2M checkpoint")

Pathways = c("[TPW] PI3Kgamma Signature", "[TBI] MAPK up genes", "[IPA] UVB Induced MAPK Signaling",
             "[IPA] ErbB2 ErbB3 Signaling", "[IPA] Estrogen Dependent Breast Cancer Signaling", "[IPA] AMPK Signaling",
             "[IPA] PTEN Signaling", "[IPA] ERK MAPK Signaling", "[HM] Estrogen response",
             "[HM] Angiogenesis", "[HM] PI3K Akt mTOR signaling", "[LM] Proliferation",
             "[IPA] Mismatch Repair in Eukaryotes", "[HM] KRAS signaling down", "[HM] E2F targets", "[HM] G2M checkpoint")

# Preparation Heatmap
Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann$IMS)),]
Clinical_data_ann = Clinical_data_ann[which(Clinical_data_ann$IMS %in% IMS_groups),]
Clinical_data_ann = Clinical_data_ann[which(Clinical_data_ann$Assigned_Ethnicity_simplified %in% c("Black", "White")),]
colnames(ES) = substring(colnames(ES), 1, 12)
ES = ES[,Clinical_data_ann$bcr_patient_barcode]

ES = ES[Pathways,]

ESz = ES
for(j in 1: nrow(ESz))  {
  ESz[j,] = (ES[j,]-mean(ES[j,]))/sd(ES[j,]) # z-score the enrichment matrix
}

Ethnicity_colors = c("White" = "#FFEC42", "Black" = "#066C3C")

IMS_colors =c("Basal" = "red", "Her2" = "pink",
              "LumA" = "darkblue", "LumB" = "lightblue",
              "Normal" = "green")
                                

#col_fun = circlize::colorRamp2(c(min(Expression.matrix.z), max(Expression.matrix.z)), c("blue", "red"))
col_fun = circlize::colorRamp2(c(-4,-0.1,0,0.1, 3), c("blue", "#EAE7EF","white","#FAD9D0", "red"))

mat =  ESz
ha_column = HeatmapAnnotation(df = data.frame(IMS = Clinical_data_ann$IMS,
                                              Ethnicity = Clinical_data_ann$Assigned_Ethnicity_simplified),
                              show_annotation_name = TRUE,
                              col = list(IMS = IMS_colors,
                                         Ethnicity = Ethnicity_colors))

## Cluster tree
col_dend = as.dendrogram(hclust(dist(t(ESz))), method = "complete")
#col_dend = color_branches(col_dend, k = 2, col = c( "purple", "orange"))
x = as.data.frame(t(ESz))
x$cluster = cutree(col_dend, k=2)[match(rownames(x), names(cutree(col_dend, k=2)))]
x$cluster[which(x$cluster == 1)] = "beneficial"
x$cluster[which(x$cluster == 2)] = "non-beneficial"

x$Ethnicity = Clinical_data_ann$Assigned_Ethnicity_simplified[match(rownames(x),
                                                                    Clinical_data_ann$bcr_patient_barcode)]

dir.create("./Analysis/Benefit_Cluster", showWarnings = FALSE)
dir.create("./Analysis/Benefit_Cluster/With_PAM50", showWarnings = FALSE)
save(x, file = paste0("./Analysis/Benefit_Cluster/With_PAM50/", paste0(IMS_groups, collapse = "_"), ".Rdata"))

# Heatmap
#png(paste0("./Figures/Heatmaps/",  Gene.set.list, "_", Gene.set,".png"), res = 600, width = 10, height = 10, units = "in")

dir.create("./Figures/Heatmaps/C_21_Heatmap_Selected_Pathways_ES/With_PAM50", showWarnings = FALSE)
png(paste0("./Figures/Heatmaps/C_21_Heatmap_Selected_Pathways_ES/With_PAM50/C_21_v2_Heatmap_Selected_Pathways_ES", 
           paste0(IMS_groups, collapse = "_"), ".png"), res = 600, width = 7, height = 5, units = "in")
plot = Heatmap(mat,
               name = "expression values", 
               show_heatmap_legend = FALSE,
               cluster_rows = TRUE,
               cluster_columns = col_dend,
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




