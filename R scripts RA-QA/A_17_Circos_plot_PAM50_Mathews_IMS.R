

## Set up environment
rm(list = ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-RA-QA")

source("../tools/ipak.function.R")
ipak(required.packages = c("circlize", "migest", "dplyr"))

## Load data
load("./Analysis/Sample_annotations.Rdata")

Clinical_data_ann = Clinical_data_ann[-which(Clinical_data_ann$Ethnicity == "Other"),]
Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann$IMS)),]
df = Clinical_data_ann[, c("IMS", "IMS_Mathews")]

plot_df = ddply(df, .(df$IMS, df$IMS_Mathews), nrow)

grid.col = c(BasalHer2 = "#F4A101", BasalLumHer2 = "#92D050", BasalMyo = "#3D64E4", Lum = "#158400", LumBasal = "#D36733",
             MyoLumA = "#6F1287", MyoLumB = "#C720F1", MyoLumHer2 = "#F377C6",
             Basal = "red", Her2 = "pink", LumA = "blue", LumB = "lightblue", Normal = "green")

dir.create("./Figures/A_17_Circos_plot", showWarnings = FALSE)
png("./Figures/A_17_Circos_plot/Circos_plot2.png", res = 600, height = 5, width = 5, units = "in")
circos.par(gap.after = c(rep(1, length(unique(df$IMS_Mathews))-1), 15, rep(1, length(unique(df$IMS))-1), 15))
chordDiagram(plot_df, order = c("BasalHer2", "BasalMyo", "BasalLumHer2", "Lum", "LumBasal", "MyoLumA", "MyoLumB", "MyoLumHer2",
                                "Basal", "Her2", "LumA", "LumB", "Normal"), grid.col = grid.col,
             annotationTrack = "grid")
#circos.track(track.index = 1, panel.fun = function(x, y) {
#  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
#              facing = "clockwise", niceFacing = TRUE, adj = c(-0.5, 0.5))
#}, bg.border = NA)
dev.off()
circos.clear()
