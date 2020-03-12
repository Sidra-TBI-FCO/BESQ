#script for ICR clustering of RNAseq data from subreads
# Setup environment
rm(list=ls())
setwd("F:/DropBox Wouter/Dropbox (TBI-Lab)/External Collaborations/NNN-BRCA/")
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA")
required.packages.BioC <- c("ConsensusClusterPlus","clue","heatmap3")
missing.packages <- required.packages.BioC[!(required.packages.BioC %in% installed.packages()[,"Package"])]
if(length(missing.packages)) {
  BiocManager::install(missing.packages)
}
library(ConsensusClusterPlus)
library(clue)
library (heatmap3)

# Set parameters
groups = "HML.ICR.Cluster" # "k4" or "HL_k4" or "k2" or "HL_k2"

# load data
load("./Data/expression matrix/BRCA_gene_RNAseq_normalized_Tissue_Type_and_whitelists_filtered.Rdata")
load("../tools/ICR_genes.RData")
#Clinical_data = read.csv("./Data/Clinical Data/TCGA_CLINICAL_DATA_CELL_2018_S1.csv", stringsAsFactors = FALSE)

# Create directories
dir.create("./Analysis/ICR data", showWarnings = FALSE)

# subset expression matix for ICR
ICR_subset_RNAseq = t(filtered.norm.RNAseqData[row.names(filtered.norm.RNAseqData) %in% ICR_genes, ])
ICR_subset_RNAseq_log2 = log(ICR_subset_RNAseq +1, 2)

# Hierarchical Clustering
source("../tools/stefanofunctions.R")
setwd("./Analysis/ICR data")

ddist = dist(ICR_subset_RNAseq_log2)

ConsensusClusterObject <- ConsensusClusterPlus(ddist,                                                                
                                               maxK = 6,                                                                              # set K
                                               pItem = 0.8,
                                               reps=5000,                                                                             # set repeats
                                               title=paste0("renormalized.Full.matrix.ICR.reps5000"),              # Output filename (no path)
                                               clusterAlg = "hc",                                                                     # clustering Algorithm : Hierarchiocal clustering
                                               innerLinkage = "ward.D2",                                                              # for color coding the clusters use tmyPal = ...
                                               finalLinkage = "complete",
                                               plot = 'pdf',                                                                          # write resut to pdf (Alt:png)
                                               writeTable = TRUE,
                                               verbose = TRUE)
outputfiles = list.files(paste0("renormalized.Full.matrix.ICR.reps5000"), full.names = TRUE)
class_files = outputfiles[grep("consensusClass", outputfiles)]

N.files = length(class_files)
table_cluster_assignment = data.frame(ICRscore = rowMeans(ICR_subset_RNAseq_log2))

for (j in 1:N.files){
  file = paste0("./", class_files[j])
  consensus_class = read.csv(file = file,header=FALSE)
  group = paste0("Group_k",j+1)
  colnames(consensus_class) = c("PatientID", group)
  rownames(consensus_class) = consensus_class$PatientID
  consensus_class$PatientID = NULL
  table_cluster_assignment[,group] = consensus_class[,group][match(rownames(table_cluster_assignment), rownames(consensus_class))]
  
  transl_table_ICR_cluster = aggregate(ICRscore~get(group),data = table_cluster_assignment, FUN=mean)
  colnames(transl_table_ICR_cluster) = c(group,"mean_ICRscore")
  transl_table_ICR_cluster = cbind(transl_table_ICR_cluster[order(transl_table_ICR_cluster$mean_ICRscore),],ICR_name=paste0("ICR",c(1:(j+1))))
  
  ICR_cluster = paste0("ICR_cluster_k",j+1)
  table_cluster_assignment[,ICR_cluster] = transl_table_ICR_cluster$ICR_name[match(table_cluster_assignment[,group],
                                                                                   transl_table_ICR_cluster[,group])]
}

#calinsky
sHc <- hclust(ddist, method = "ward.D2")
aCalinsky <- calinsky(sHc,gMax=10)
pdf(file = paste0("./renormalized.Full.matrix.ICR.reps5000/ICR_cluster_assignment_k2-6.Calinsky.pdf"), width = 16, height = 6)
plot(aCalinsky, type = "l", col = "grey", main = "Calinsky & Harabasz curve", xlab = "# of groups")
text(1:length(aCalinsky), aCalinsky, paste(1:length(aCalinsky)))
dev.off()
optimal.calinsky = which(aCalinsky == max(aCalinsky[3:5]))

clustering = table_cluster_assignment

#save data
outputname = paste0("./TCGA_BRCA_table_cluster_assignment.RData")
save(clustering,optimal.calinsky, file = outputname) 

load("./Analysis/ICR data/TCGA_BRCA_table_cluster_assignment.RData")

# HML assignment
clustering$HML.ICR.Cluster = NA
clustering$HML.ICR.Cluster[which(clustering$ICR_cluster_k4 == "ICR1")] = "ICR Low"
clustering$HML.ICR.Cluster[which(clustering$ICR_cluster_k4 == "ICR2")] = "ICR Medium"
clustering$HML.ICR.Cluster[which(clustering$ICR_cluster_k4 == "ICR3")] = "ICR Medium"
clustering$HML.ICR.Cluster[which(clustering$ICR_cluster_k4 == "ICR4")] = "ICR High"

save(clustering,optimal.calinsky, file = outputname) 

setwd("../../")
load("./Analysis/ICR data/TCGA_BRCA_table_cluster_assignment.RData")

if(combine_ICR_Medium_low = "combine_ICR_Medium_low"){
  clustering$HML.ICR.Cluster[which(clustering$HML.ICR.Cluster %in% c("ICR Low", "ICR Medium"))] = "ICR Medium-Low"
}


#ICR heatmap
dir.create("./Figures/Heatmaps", showWarnings = FALSE)
dir.create("./Figures/Heatmaps/ICR_heatmaps", showWarnings = FALSE)

clustering <- clustering[order(clustering$ICRscore),]
expression.matrix.subset.ICR <- filtered.norm.RNAseqData[ICR_genes,rownames(clustering)]
colors = clustering[,c("ICRscore", "HML.ICR.Cluster")]

if(groups == "HML.ICR.Cluster"){
  colors$HML.ICR.Cluster[colors$HML.ICR.Cluster=="ICR Low"]<-"blue"
  colors$HML.ICR.Cluster[colors$HML.ICR.Cluster=="ICR Medium"]<-"green"
  colors$HML.ICR.Cluster[colors$HML.ICR.Cluster=="ICR High"]<-"red"
  colors$HML.ICR.Cluster[colors$HML.ICR.Cluster=="ICR Medium-Low"]<-"purple"
}

#colors$Grade[colors$Grade==2]<-"pink"
#colors$Grade[colors$Grade==3]<-"purple"
#colnames(colors) <- c("Grade","Type","HL.Cluster")
#colors = colors[order(colors$HL.Cluster),]
colors <- colors[,2,drop=FALSE]
colors = as.matrix(colors)

# reoder expression matrix based on colors
ICR_subset_RNAseq_log2 = t(ICR_subset_RNAseq_log2)
expression.matrix.subset.ICR = ICR_subset_RNAseq_log2[rev(ICR_genes),rownames(colors)]

png(filename = paste0("./Figures/Heatmaps/ICR_heatmaps/ICR_heatmapv3_", groups, ".png"), res = 600,
    width = 4, height = 4.5, units = "in")
heatmap3(expression.matrix.subset.ICR,
         main="TCGA BRCA ICR RNASeq",
         ColSideColors=colors,
         Colv= NA , 
         Rowv = NA,
         #as.dendrogram(sHc),
         #col=bluered(75) ,
         #labRow = NA,
         #labCol = samples,
         scale='row',
         margins = c(12, 7),
         labCol = NA)
#legend("topright",legend = c("ICR Low", "ICR Medium 1", "ICR Medium 2", "ICR High","","TNBC","Non-TNBC","","Grade 2","Grade 3"),
       #col = c("blue", "green", "yellow" ,"red","white","orange","lightblue","white","pink","purple"),lty= 1,lwd = 5,cex = 0.7)

dev.off()
