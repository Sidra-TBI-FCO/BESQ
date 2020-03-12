#######
#
# Heatmap of PAM50 genes by Mathews et al paper
#
#######


# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA")
source("../tools/ipak.function.R")
required.packages <- c("ComplexHeatmap", "impute", "ctc", "amap")
ipak(required.packages)

#Set parameters

# Load data
load("../tools/Selected.pathways.3.4.RData")
load("./Analysis/Sample_annotations.Rdata")
load("./Data/expression matrix/BRCA_gene_RNAseq_normalized_Tissue_Type_and_whitelists_filtered.Rdata")

# Preparation Heatmap
Clinical_data_ann = Clinical_data_ann[which(Clinical_data_ann$Assigned_Ethnicity_simplified %in% c("White", "Black")),]
Clinical_data_ann = Clinical_data_ann[which(Clinical_data_ann$IMS_Mathews %in% c("BasalMyo")),]
Clinical_data_ann = Clinical_data_ann[order(Clinical_data_ann$IMS),]
Clinical_data_ann$Ethnicity_ICR = paste(Clinical_data_ann$Assigned_Ethnicity_simplified, Clinical_data_ann$HML.ICR.Cluster)

RNASeq.QN.LOG2 = log(filtered.norm.RNAseqData + 1, 2)

selected_genes = unname(unlist(Selected.pathways$`[IPA] AMPK Signaling`))
selected_genes = selected_genes[which(selected_genes %in% rownames(RNASeq.QN.LOG2))]
selected_genes = selected_genes[-which(selected_genes == "INS")]

colnames(RNASeq.QN.LOG2) = substring(colnames(RNASeq.QN.LOG2), 1, 12)

Expression.matrix = RNASeq.QN.LOG2[selected_genes,Clinical_data_ann$bcr_patient_barcode]

Expression.matrix.z = Expression.matrix
for(j in 1: nrow(Expression.matrix.z))  {
  Expression.matrix.z[j,] = (Expression.matrix[j,]-mean(Expression.matrix[j,]))/sd(Expression.matrix[j,]) # z-score the enrichment matrix
}

#col_fun = circlize::colorRamp2(c(min(Expression.matrix.z), max(Expression.matrix.z)), c("blue", "red"))
col_fun = circlize::colorRamp2(c(-4,-0.1,0,0.1, 7), c("blue", "#EAE7EF","white","#FAD9D0", "red"))

Ethnicity_Black_ICR_colors = c("Black ICR High" = "red", "Black ICR Low" = "blue", "Black ICR Medium" = "green",
                               "White ICR High" = "grey", "White ICR Low" = "grey", "White ICR Medium" = "grey")

Ethnicity_White_ICR_colors = c("Black ICR High" = "grey", "Black ICR Low" = "grey", "Black ICR Medium" = "grey",
                               "White ICR High" = "red", "White ICR Low" = "blue", "White ICR Medium" = "green")

Ethnicity_colors = c("White" = "#FFEC42", "Black" = "#066C3C")

ha_column = HeatmapAnnotation(df = data.frame(Ethnicity = Clinical_data_ann$Assigned_Ethnicity_simplified,
                                              Ethnicity_ICR_Black = Clinical_data_ann$Ethnicity_ICR,
                                              Ethnicity_ICR_White = Clinical_data_ann$Ethnicity_ICR),
                              show_annotation_name = TRUE,
                              col = list(Ethnicity = Ethnicity_colors,
                                         Ethnicity_ICR_Black = Ethnicity_Black_ICR_colors,
                                         Ethnicity_ICR_White = Ethnicity_White_ICR_colors))
mat = Expression.matrix.z[, which(substring(colnames(Expression.matrix.z), 1, 12) %in% 
                                    Clinical_data_ann$bcr_patient_barcode)]

plot = Heatmap(mat,
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.rect(x = x, y = y, width = width, height = height, gp = gpar(fill =col_fun(mat[i, j]),lty=0))},
               name = "expression values", 
               show_heatmap_legend = FALSE,
               cluster_rows = TRUE,
               cluster_columns = TRUE,
               row_title_gp = gpar(fontsize = 0.1),
               #column_names_gp = gpar(fontsize = 10),
               row_names_gp = gpar(fontsize = 10),
               top_annotation = ha_column,
               #column_title = paste0("Expression PAM50 genes"),
               show_column_names = FALSE,
               show_row_names = TRUE)
# Heatmap
#png(paste0("./Figures/Heatmaps/",  Gene.set.list, "_", Gene.set,".png"), res = 600, width = 10, height = 10, units = "in")

dir.create("./Figures/Heatmaps/AMPK_Heatmaps/", showWarnings = FALSE)
png(paste0("./Figures/Heatmaps/AMPK_Heatmaps/AMPK_Heatmap_v1.png"), res = 600, width = 10, height = 10, units = "in")

draw(plot, show_annotation_legend = FALSE)

dev.off()




