
# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-RA-QA")                                                                   # Setwd to location were output files have to be saved.

source("../tools/ipak.function.R")

ipak(c("ggplot2", "ggpubr"))

required.bioconductor.packages = c("GSVA","heatmap3", "gclus")                                                                   
ibiopak(required.bioconductor.packages)

# Set parameters
gene_list = "Bindea_ORIG"

# Load data
load("../tools/ICR_genes.RData")
load("./Data/expression matrix/CTA_TNBC_BRCA_Qatar_RNASeq.tumor.tissue.filtered.EDAseq.QN.Rdata")
load("../tools/Selected.pathways.3.4.RData")
load("../tools/immune.gene.lists.v3.Rdata")

Selected.pathways[[55]] = ICR_genes
names(Selected.pathways)[55] = "ICR genes"

Expression.data = RNASeq.QN.LOG2

available_genes = rownames(Expression.data)
Gene.list = get(gene_list)
unavailable_genes_RNAseq = unlist(Gene.list)[-which(unlist(Gene.list) %in% rownames(Gene.list))]

## ssGSEA
ES = gsva(Expression.data,Gene.list,method="ssgsea")

dir.create("./Analysis/ssGSEA_scores", showWarnings = FALSE)
save(file = paste0("./Analysis/ssGSEA_scores/", gene_list, "_ES_scores.Rdata"), ES)
