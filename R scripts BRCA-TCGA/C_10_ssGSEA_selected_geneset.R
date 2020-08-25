
# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA")                                                                   # Setwd to location were output files have to be saved.

source("../tools/ipak.function.R")

ipak(c("ggplot2", "ggpubr"))

required.bioconductor.packages = c("GSVA","heatmap3", "gclus")                                                                   
ibiopak(required.bioconductor.packages)

# Set parameters
gene_list = "Selected.pathways"

# Load data
load("../tools/ICR_genes.RData")
load("./Data/expression matrix/BRCA_gene_RNAseq_normalized_Tissue_Type_and_whitelists_filtered.Rdata")
load("../tools/Selected.pathways.3.4.RData")
load("../tools/immune.gene.lists.v3.Rdata")

Selected.pathways[[55]] = ICR_genes
names(Selected.pathways)[55] = "ICR genes"

Selected.pathways[[56]] = c("ERBB2", "ERBB3", "ERBB4")
names(Selected.pathways)[56] = "ERBB2, ERBB3, and ERBB4"

RNAseq_log2 = log(filtered.norm.RNAseqData +1, 2)
Expression.data = RNAseq_log2

available_genes = rownames(Expression.data)
Gene.list = get(gene_list)
unavailable_genes_RNAseq = unlist(Gene.list)[-which(unlist(Gene.list) %in% rownames(Gene.list))]

## ssGSEA
ES = gsva(Expression.data,Gene.list,method="ssgsea")

dir.create("./Analysis/ssGSEA_scores", showWarnings = FALSE)
save(file = paste0("./Analysis/ssGSEA_scores/v2_", gene_list, "_ES_scores.Rdata"), ES)
