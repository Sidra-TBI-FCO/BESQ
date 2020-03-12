#######
#
# Correlation matrix between genes of interest
# PAM50 obtained from: https://www.biostars.org/p/210106/
#
#######


# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA")
#setwd("~/Dropbox (TBI-Lab)/External Collaborations/NNN-BRCA/RNAseq-TCGA_BRCA")
source("../tools/ipak.function.R")
required.packages <- c()
ipak(required.packages)  

#source("./Analysis QBRI and TCGA/R tools/PAM50_Parker/")

#Set parameters

# Load data
load("./Data/expression matrix/BRCA_gene_RNAseq_normalized_Tissue_Type_and_whitelists_filtered.Rdata")

# PAM50 classification
## RNASeq Data from the EDASeq protocol after log2 transformation
RNASeq.QN.LOG2 = log(filtered.norm.RNAseqData +1, 2)
dir.create("./Analysis/PAM50", showWarnings = FALSE)

Basal.genes=c("PTTG1", "CDC20", "ORC6", "KIF2C", "UBE2C", "MELK", "BIRC5",
             "NUF2", "CEP55", "EXO1", "CENPF", "NDC80", "TYMS", "UBE2T",
             "ANLN", "CCNB1", "RRM2", "MKI67")
Myo1.genes=c("MIA", "FOXC1", "ACTR3B", "CCNE1", "PHGDH")
Myo2.genes=c("SFRP1", "MYC", "CDH3", "KRT5", "KRT17", "KRT14", "EGFR")
Luminal.genes=c("CDC6", "ESR1", "BCL2", "FOXA1", "CXXC5", "MLPH", "MAPT",
                "SLC39A6", "NAT1", "MDM2", "PGR", "MMP11", "BLVRA")
Her2.genes=c("FGFR4", "GRB7", "ERBB2")

PAM50_genes = c(Basal.genes,Myo1.genes,Myo2.genes,Luminal.genes,Her2.genes)

PAM50_table = data.frame(Gene = PAM50_genes, Category = c(rep("basal", length(Basal.genes)),
                                                          rep("myo1", length(Myo1.genes)),
                                                          rep("myo2", length(Myo2.genes)),
                                                          rep("luminal", length(Luminal.genes)),
                                                          rep("Her2", length(Her2.genes))))
PAM50_table$Category = factor(PAM50_table$Category, levels = c("basal", "myo1", "myo2", "luminal", "Her2"))

Subtypes = c("Basal","Her2","LumA","LumB","Normal")

save(PAM50_genes, PAM50_table, file = "../tools/PAM50_genes.RData")

#Set parameters
IMS_colors = list(IMS_PAM50 = c("Basal" = "red", "Her2" = "pink", 
                                "LumA" = "darkblue", "LumB" = "lightblue",
                                "Normal" = "green"))

PAM50.matrix=t(RNASeq.QN.LOG2[PAM50_genes,])
write.csv(PAM50.matrix,file = "./Analysis/PAM50/PAM50_matrix_Log2.csv")

#Wouts Method
Mean.expressions = data.frame(Basal.sig = colMeans(RNASeq.QN.LOG2[Basal.genes,]),
                              Myo1.sig = colMeans(RNASeq.QN.LOG2[Myo1.genes,]),
                              Myo2.sig = colMeans(RNASeq.QN.LOG2[Myo2.genes,]),
                              Luminal.sig = colMeans(RNASeq.QN.LOG2[Luminal.genes,]),
                              Her2.sig = colMeans(RNASeq.QN.LOG2[Her2.genes,]))

computation.matrix <- read.csv("../tools/PAM_50_Mathews/signature_table.csv")
# name a b c d e
#   1     BasalMyo + + + - -
#   2    BasalHer2 + B B - +
#   3 BasalLumHer2 + B B + +
#   4          Lum - - - + -
#   5     LumBasal + - - + -
#   6   MyoLumHer2 - B + B +
#   7      MyoLumB - + + B B
#   8      MyoLumA - - + B B

Mathews.scores = data.frame (BasalMyo = Mean.expressions$Basal.sig + Mean.expressions$Myo1.sig + Mean.expressions$Myo2.sig - Mean.expressions$Luminal.sig - Mean.expressions$Her2.sig,
                             BasalHer2 = Mean.expressions$Basal.sig - Mean.expressions$Luminal.sig + Mean.expressions$Her2.sig,
                             BasalLumHer2 = Mean.expressions$Basal.sig + Mean.expressions$Luminal.sig - Mean.expressions$Her2.sig,
                             Lum = - Mean.expressions$Basal.sig - Mean.expressions$Myo1.sig - Mean.expressions$Myo2.sig + Mean.expressions$Luminal.sig - Mean.expressions$Her2.sig)


#Mathews Method
data_file_metabric <- "../tools/PAM_50_Mathews/data_pam46_metabric.csv"
data_file_TCGA <- "./Analysis/PAM50/PAM50_matrix_Log2.csv"
feature_group_file <- "../tools/PAM_50_Mathews/feature_groups.csv"
signature_table_file <- "../tools/PAM_50_Mathews/signature_table.csv"

zscore_normalize <- function(sample_data) {
  for(i in 1:ncol(sample_data)) {
    variance_i <- var(sample_data[,i])
    mean_i <- mean(sample_data[,i])
    sample_data[,i] <- (sample_data[,i]-mean_i) / sqrt(variance_i)
  }
  return(sample_data)
}

adjust_value_by_table_entry <- function(value, symbol) {
  if(symbol == "+") {
    return(value)
  }
  if(symbol == "-") {
    return(-1*value)
  }
  if(symbol == "B") {
    return(0)
  }
  return(NA)
}

distance_to_signature <- function(sample, signature, feature_groups) {
  group_means <- lapply(feature_groups, function(x) {mean(sample[x[,1],])} )
  
  tally <- 0
  for(i in 1:length(signature)) {
    tally <- tally + adjust_value_by_table_entry(group_means[[i]], signature[i])
  }
  # Exponential is just to make something like a 'distance' rather than a score (this function is monotone, and does not affect the classification)
  return(exp(-1*tally))
}

get_classification <- function(data_file, feature_group_file, signature_table_file) {
  feature_groups <- read.csv(feature_group_file, header=FALSE, stringsAsFactors=FALSE)
  feature_groups <- split(feature_groups, feature_groups[,2])
  signature_table <- read.csv(signature_table_file, stringsAsFactors=FALSE)
  rownames(signature_table) <- signature_table$name
  signature_table$name <- c()
  signature_table <- t(signature_table)
  
  sample_data <- read.csv(data_file, stringsAsFactors=FALSE)
  rownames(sample_data) <- sample_data[,1]
  sample_data[,1] <- c()
  sample_data <- zscore_normalize(sample_data)
  
  # Finds a 'distance' associated with each signature and each sample
  number_signatures <- dim(signature_table)[2]
  distances <- c()
  for(i in 1:number_signatures) {
    signature <- signature_table[,i]
    distances_to_signature <- apply(sample_data, 1, function(sample){distance_to_signature(data.frame(sample), signature, feature_groups)} )
    distances <- cbind(distances, distances_to_signature)
  }
  distances <- data.frame(distances)
  colnames(distances) <- colnames(signature_table)
  
  # For a given sample chooses the class with the smallest 'distance'
  classes <- data.frame(apply(distances, 1, function(x) {colnames(distances)[which.min(x)]} ))    
  colnames(classes) <- c("classes")
  return(classes)
}

classes <- get_classification(data_file_TCGA, feature_group_file, signature_table_file)
write.csv(classes, "./Analysis/PAM50/Mathews_classification.csv", quote=FALSE)
