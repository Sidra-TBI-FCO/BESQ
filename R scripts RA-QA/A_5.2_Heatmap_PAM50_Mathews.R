#######
#
# Heatmap of PAM50 genes by Mathews et al paper
#
#######


# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-RA-QA")
source("../tools/ipak.function.R")
required.packages <- c("ComplexHeatmap","amap")
ipak(required.packages)  

# Load data
load("../tools/PAM50_genes.RData")
load("./Analysis/Sample_annotations.Rdata")
load("./Data/expression matrix/CTA_TNBC_BRCA_Qatar_RNASeq.tumor.tissue.filtered.EDAseq.QN.Rdata")


# Preparation Heatmap
Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann$IMS_Mathews)),]
Clinical_data_ann = Clinical_data_ann[order(Clinical_data_ann$IMS_Mathews),]

Expression.matrix = RNASeq.QN.LOG2[PAM50_genes,]

Expression.matrix.z = Expression.matrix
for(j in 1: nrow(Expression.matrix.z))  {
  Expression.matrix.z[j,] = (Expression.matrix[j,]-mean(Expression.matrix[j,]))/sd(Expression.matrix[j,]) # z-score the enrichment matrix
}

IMS_Mathews_colors = c("BasalHer2" = "#F4A101", 
                       "BasalMyo" = "#3D64E4",
                       "BasalLumHer2" = "#92D050",
                       "MyoLumB" = "#C720F1",
                       "Lum" = "#158400","LumBasal" = "#D36733",
                       "MyoLumA" = "#6F1287", "MyoLumHer2" = "#F377C6")

IMS_colors = c("Basal" = "red", "Her2" = "pink", 
               "LumA" = "darkblue", "LumB" = "lightblue",
               "Normal" = "green")

#col_fun = circlize::colorRamp2(c(min(Expression.matrix.z), max(Expression.matrix.z)), c("blue", "red"))
col_fun = circlize::colorRamp2(c(-3,-0.1,0,0.1, 2.5), c("blue", "#EAE7EF","white","#FAD9D0", "red"))
Subtypes = unique(Clinical_data_ann$IMS_Mathews)

for (i in 1:length(Subtypes)){
  subtype = Subtypes[i]
  assign(paste0("mat", i), Expression.matrix.z[, which(substring(colnames(Expression.matrix.z), 1, 12) %in% 
                                                         Clinical_data_ann$Sample_ID[which(Clinical_data_ann$IMS_Mathews == subtype)]),drop=FALSE])
  assign(paste0("ha_column", i), HeatmapAnnotation(df = data.frame(IMS_Mathews = Clinical_data_ann$IMS_Mathews[which(Clinical_data_ann$IMS_Mathews == subtype)],
                                                                   IMS_PAM50 = Clinical_data_ann$IMS[which(Clinical_data_ann$IMS_Mathews == subtype)]),
                                                   show_annotation_name = FALSE,
                                                   col = list(IMS_Mathews = IMS_Mathews_colors,
                                                              IMS_PAM50 = IMS_colors)))
}

# Heatmap
#png(paste0("./Figures/Heatmaps/",  Gene.set.list, "_", Gene.set,".png"), res = 600, width = 10, height = 10, units = "in")

dir.create("./Figures/Heatmaps/PAM50_Heatmaps_Mathews/", showWarnings = FALSE)
png(paste0("./Figures/Heatmaps/PAM50_Heatmaps_Mathews/PAM50_Heatmap_v4.png"), res = 600, width = 6, height = 7.5, units = "in")
plot = Heatmap(mat1,
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.rect(x = x, y = y, width = width, height = height, gp = gpar(fill =col_fun(mat1[i, j]),lty=0))},
               name = "expression values", 
               show_heatmap_legend = FALSE,
               cluster_rows = FALSE,
               #cluster_columns = FALSE,
               row_title_gp = gpar(fontsize = 0.1),
               #column_names_gp = gpar(fontsize = 10),
               row_names_gp = gpar(fontsize = 10),
               top_annotation = ha_column1,
               #column_title = paste0("Expression PAM50 genes"),
               show_column_names = FALSE,
               show_row_names = FALSE,split = PAM50_table$Category
) +
  Heatmap(mat2,
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.rect(x = x, y = y, width = width, height = height, gp = gpar(fill =col_fun(mat2[i, j]),lty=0))},
          top_annotation = ha_column2, show_column_names = FALSE,
          show_row_names = FALSE, show_heatmap_legend = FALSE, cluster_rows = FALSE) +
  Heatmap(mat3,
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.rect(x = x, y = y, width = width, height = height, gp = gpar(fill =col_fun(mat3[i, j]),lty=0))},
          top_annotation = ha_column3, show_column_names = FALSE,
          show_row_names = FALSE, show_heatmap_legend = FALSE, cluster_rows = FALSE) +
  Heatmap(mat4,
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.rect(x = x, y = y, width = width, height = height, gp = gpar(fill =col_fun(mat4[i, j]),lty=0))},
          top_annotation = ha_column4, show_column_names = FALSE,
          show_row_names = FALSE, show_heatmap_legend = FALSE, cluster_rows = FALSE) +
  Heatmap(mat5,
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.rect(x = x, y = y, width = width, height = height, gp = gpar(fill =col_fun(mat5[i, j]),lty=0))},
          top_annotation = ha_column5, show_column_names = FALSE,
          show_row_names = FALSE, show_heatmap_legend = FALSE, cluster_rows = FALSE) +
  Heatmap(mat6,
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.rect(x = x, y = y, width = width, height = height, gp = gpar(fill =col_fun(mat6[i, j]),lty=0))},
          top_annotation = ha_column6, show_column_names = FALSE,
          show_row_names = FALSE, show_heatmap_legend = FALSE, cluster_rows = FALSE) +
  Heatmap(mat7,
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.rect(x = x, y = y, width = width, height = height, gp = gpar(fill =col_fun(mat7[i, j]),lty=0))},
          top_annotation = ha_column7, show_column_names = FALSE,
          show_row_names = FALSE, show_heatmap_legend = FALSE, cluster_rows = FALSE) +
  Heatmap(mat8,
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.rect(x = x, y = y, width = width, height = height, gp = gpar(fill =col_fun(mat8[i, j]),lty=0))},
          top_annotation = ha_column8, show_column_names = FALSE,
          show_row_names = TRUE, show_heatmap_legend = FALSE, cluster_rows = TRUE)

draw(plot, show_annotation_legend = FALSE, gap = unit(1, "mm"))
dev.off()




