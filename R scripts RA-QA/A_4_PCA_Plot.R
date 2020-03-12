#######
#
# PCA plot of samples
#
#######


# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-RA-QA")
source("../tools/ipak.function.R")
required.packages <- c("corrplot", "stringr")
ipak(required.packages)  
required.packages = c("ggrepel", "ggplot2")
ipak(required.packages)

# Set parameters
selected_genes = "ICR"
Type = "Normalized" # "Raw" or "Normalized"

# Create directories
dir.create("./Figures/PCA Plots", showWarnings = FALSE)

# load data
load("./Data/expression matrix/RNASeq.filtered.tissue.raw.Rdata")
load("./Data/expression matrix/CTA_TNBC_BRCA_Qatar_RNASeq.tumor.tissue.filtered.EDAseq.QN.Rdata")
load("./Analysis/ICR data/CTA_TNBC_BRCA_Qatar_RNASeq.tumor.Clustering.RData")

#PCA
if(Type == "Raw"){
  Expression.matrix = t(expression.filtered.tissue)
}
if(Type == "Normalized"){
  Expression.matrix = t(RNASeq.NORM.quantiles)
}


#Delete columns (genes) with a value of zero for all samples (prcomp cannot create a pca object if there are constant variables)
Expression.matrix <- Expression.matrix[,apply(Expression.matrix, 2, var, na.rm=TRUE) != 0]

pca = prcomp(Expression.matrix, center= TRUE, scale=TRUE)

scores = data.frame(pca$x[,1:3])

scores = scores[match(rownames(Expression.matrix), rownames(scores)),]
scores$ICR = clustering$HL.ICR.Cluster[match(rownames(Expression.matrix), clustering$Sample_ID)]

scores$ICR[scores$ICR == "ICR-High"] = "red"
scores$ICR[scores$ICR == "ICR-Low"] = "blue"

#c(2,3,4,7,8,9,12,13,14)
p.variance.explained = round(100*(pca$sdev^2 / sum(pca$sdev^2)),1)
## 3Dplot
png(paste0("./Figures/PCA Plots/", Type, "_All_samples_by_ICR_complete_matrix.png"),res=600,height=6,width=7,unit="in")
pcaplot <-scatterplot3d(scores[,c(1:3)],
                        pch = par("pch"),
                        cex.symbols=3,
                        color=scores$ICR,
                        type = "h",
                        angle = 20,
                        main="3-D PCA results from different samples",
                        xlab=paste0("PC1 (",p.variance.explained[1],"%)"),
                        ylab=paste0("PC2 (",p.variance.explained[2],"%)"),
                        zlab=paste0("PC3 (",p.variance.explained[3],"%)")
)
plot.coords <- pcaplot$xyz.convert(scores[,c(1:3)])   # convert 3D coords to 2D projection
text(plot.coords$x, plot.coords$y,labels=rownames(scores),cex=.9, pos=4)

#legend("bottomright", 
# legend = c("C","G", "S"),
# col = c("#F45400", "#0000AF", "#FFDE02"),
# title = "Region",
# pch = 16,
# cex = 1.3,
# inset = 0, xpd = TRUE#, horiz = TRUE
#)

dev.off()
summary(pca)

scores$ICR[scores$ICR == "red"] = "ICR-High"
scores$ICR[scores$ICR == "blue"] = "ICR-Low"

png(paste0("./Figures/PCA Plots/", Type,"_2D_All_samples_by_ICR_complete_matrix.png"),res=600,height=6,width=8,unit="in")
ggplot(scores, aes(x=PC1, y=PC2,label=rownames(scores), color = ICR)) +
  geom_point(size=2, aes(col=ICR))+
  geom_text_repel() +
  #scale_colour_brewer(palette = "Set2")+
  scale_color_manual(name="ICR", values = c("red","blue"))+
  labs(x=paste0("PC1 (",p.variance.explained[1],"%)"),y=paste0("PC2 (",p.variance.explained[2],"%)"))+ 
  #scale_color_discrete(name = "ICR")+
  stat_ellipse()+  # PC1&PC2 colour by group  +
  theme_light(base_size=20) 
dev.off()

