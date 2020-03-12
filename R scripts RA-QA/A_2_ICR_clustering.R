#script for ICR clustering of RNAseq data from subreads
# Setup environment
rm(list=ls())
setwd("F:/DropBox Wouter/Dropbox (TBI-Lab)/External Collaborations/NNN-BRCA/")
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-RA-QA")
required.packages.BioC <- c("ConsensusClusterPlus","clue","heatmap3")
missing.packages <- required.packages.BioC[!(required.packages.BioC %in% installed.packages()[,"Package"])]
if(length(missing.packages)) {
  BiocManager::install(missing.packages)
}
library(ConsensusClusterPlus)
library(clue)
library (heatmap3)

# Set parameters
groups = "HL_k4" # "k4" or "HL_k4" or "k2" or "HL_k2"

# load data
load("./Data/expression matrix/CTA_TNBC_BRCA_Qatar_RNASeq.tumor.tissue.filtered.EDAseq.QN.Rdata")

# annotation data
annotation = read.csv("./Data/Clinical Data/Julie cases_clinical data_v1.csv", stringsAsFactors = FALSE)
annotation$Coding <- gsub(annotation$Coding,pattern = "-",replacement = ".")
annotation$Coding <- paste0("S.",annotation$Coding)
annotation <- annotation[annotation$Coding %in% colnames(RNASeq.NORM.quantiles),]
rownames(annotation) <- annotation$Coding
annotation$X[annotation$X ==""] <- "Non-TNBC"
RNASeq.NORM.quantiles <- RNASeq.NORM.quantiles[,annotation$Coding]
RNASeq.QN.LOG2 <- log(RNASeq.NORM.quantiles+1,2)
save(RNASeq.NORM.quantiles,RNASeq.QN.LOG2,annotation,file = "./Data/expression matrix/CTA_TNBC_BRCA_Qatar_RNASeq.tumor.tissue.annotated.EDAseq.QN.Rdata")

# subset expression matix for ICR
Gene.database <- read.csv("../tools/Gene_selection_v2.7.txt",stringsAsFactors = FALSE)
ICR.genes = Gene.database$HUGO.GENE.SYM[Gene.database$DBGS3 == 1]
ICR.genes.available = ICR.genes[which(ICR.genes %in% rownames(RNASeq.QN.LOG2))]
ICR.genes.unavailable = ICR.genes[-which(ICR.genes %in% rownames(RNASeq.QN.LOG2))]
expression.matrix.subset.ICR = RNASeq.QN.LOG2[ICR.genes.available,]

# Hierarchical Clustering
source("../tools/stefanofunctions.R")
sHc <- hclust(ddist <- dist(t(expression.matrix.subset.ICR)), method = "ward.D2")
aCalinsky <- calinsky(sHc,gMax=10)
dev.new()
plot(aCalinsky, type = "l", col = "grey", main = "Calinsky & Harabasz curve", xlab = "# of groups")
text(1:length(aCalinsky), aCalinsky, paste(1:length(aCalinsky)))
setwd("./Analysis/ICR data")
ConsensusClusterObject <- ConsensusClusterPlus(ddist,                                                                
                                               maxK = 5,                                                                          # set K
                                               pItem = 0.8,
                                               reps=5000,                                                                         # set repeats
                                               title=paste0("renormalized.Full.matrix.ICR.reps5000"),   # Output filename (no path)
                                               clusterAlg = "hc",                                                   # clustering Algorithm : Hierarchiocal clustering
                                               innerLinkage = "ward.D2",                                            # for color coding the clusters use tmyPal = ...
                                               finalLinkage = "complete",
                                               plot = 'pdf',                                                        # write resut to pdf (Alt:png)
                                               writeTable = TRUE,
                                               verbose = TRUE)

setwd("../../")
if(groups %in% c("k4", "HL_k4")){
  clustering  <- read.csv (paste0("./Analysis/ICR data/renormalized.Full.matrix.ICR.reps5000/renormalized.Full.matrix.ICR.reps5000.k=4.consensusClass.csv"),stringsAsFactors = FALSE,header = FALSE)
}
if(groups %in% c("k2", "HL_k2")){
  clustering  <- read.csv (paste0("./Analysis/ICR data/renormalized.Full.matrix.ICR.reps5000/renormalized.Full.matrix.ICR.reps5000.k=2.consensusClass.csv"),stringsAsFactors = FALSE,header = FALSE)
}

colnames(clustering) <- c("Sample_ID","Cluster")

#ICR score
ICR.scores <- as.data.frame(colMeans(expression.matrix.subset.ICR))
colnames(ICR.scores) <- "ICR.score"
clustering$ICR.score <- ICR.scores$ICR.score[match(clustering$Sample_ID,rownames(ICR.scores))]
mean.ICR.score <- aggregate(clustering$ICR.score,by = list(clustering$Cluster),FUN= mean)
mean.ICR.score <- mean.ICR.score[order(x = mean.ICR.score$x),]
if(groups %in% c("k4", "HL_k4")){
  mean.ICR.score$Cluster = paste0("ICR",c(1:4))
}
if(groups == "k4"){
  mean.ICR.score$Cluster.HL <- paste0("ICR",c("-1","-2","-3","-4"))
}
if(groups == "HL_k4"){
  mean.ICR.score$Cluster.HL <- paste0("ICR",c("-Low","-Low","-Low","-High"))
}
if(groups %in% c("k2", "HL_k2")){
  mean.ICR.score$Cluster = paste0("ICR", c(1, 2))
}
if(groups == "HL_k2"){
  mean.ICR.score$Cluster.HL <- paste0("ICR",c("-Low","-High"))
}

# Load data
load("./Analysis/ICR data/CTA_TNBC_BRCA_Qatar_RNASeq.tumor.Clustering.RData")
load("../tools/ICR_genes.RData")

clustering$HL.ICR.Cluster <- mean.ICR.score$Cluster.HL[match(clustering$Cluster,mean.ICR.score$Group.1)] 


#ICR heatmap
dir.create("./Figures/Heatmaps/ICR_heatmaps", showWarnings = FALSE)

clustering <- clustering[order(clustering$ICR.score),]
expression.matrix.subset.ICR <- expression.matrix.subset.ICR[,clustering$Sample_ID]
colors <- annotation[clustering$Sample_ID,c("Grade","X")]
colors$HL.Cluster <- clustering$HL.ICR.Cluster[match(rownames(colors),clustering$Sample_ID)]
colors$X[colors$X=="TNBC"]<-"orange"
colors$X[colors$X=="Non-TNBC"]<-"lightblue"



if(groups == "k4"){
  colors$HL.Cluster[colors$HL.Cluster=="ICR-1"]<-"blue"
  colors$HL.Cluster[colors$HL.Cluster=="ICR-2"]<-"green"
  colors$HL.Cluster[colors$HL.Cluster=="ICR-3"]<-"yellow"
  colors$HL.Cluster[colors$HL.Cluster=="ICR-4"]<-"red"
}
if(groups %in% c("HL_k4", "HL_k2")){
  colors$HL.Cluster[colors$HL.Cluster=="ICR-Low"]<-"blue"
  colors$HL.Cluster[colors$HL.Cluster=="ICR-High"]<-"red"
}

colors$Grade[colors$Grade==2]<-"pink"
colors$Grade[colors$Grade==3]<-"purple"
colnames(colors) <- c("Grade","Type","HL.Cluster")
colors = colors[order(colors$HL.Cluster),]
colors$Grade = NULL
colors$Type = NULL
colors <- as.matrix(colors)


test = rev(ICR_genes)
# reoder expression matrix based on colors
expression.matrix.subset.ICR = expression.matrix.subset.ICR[test,rownames(colors)]

png(filename = paste0("./Figures/Heatmaps/ICR_heatmaps/ICR_heatmapv2_", groups, ".png"), res = 600,
    width = 6, height = 6, units = "in")
heatmap3(expression.matrix.subset.ICR,
         main="RA-QA ICR RNASEQ",
         ColSideColors=colors,
         Colv= NA,
         Rowv = NA,
         #as.dendrogram(sHc),
         #col=bluered(75) ,
         labCol = NA,
         #labCol = samples,
         scale='row',
         margins = c(12, 7))
#legend("topright",legend = c("ICR Low", "ICR Medium 1", "ICR Medium 2", "ICR High","","TNBC","Non-TNBC","","Grade 2","Grade 3"),
       #col = c("blue", "green", "yellow" ,"red","white","orange","lightblue","white","pink","purple"),lty= 1,lwd = 5,cex = 0.7)

dev.off()
save(expression.matrix.subset.ICR,colors,clustering,file = "./RNAseq Data/ICR data/CTA_TNBC_BRCA_Qatar_RNASeq.tumor.Clustering.RData")
