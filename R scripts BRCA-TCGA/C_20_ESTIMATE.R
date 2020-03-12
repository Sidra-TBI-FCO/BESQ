#################################################################
###
### 
#################################################################

# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA") 

source("~/Dropbox (TBI-Lab)/DB-LAB/bioinformatics tools/R scripts/ipak.function.R")

required.packages <- c("corrplot", "stringr")
ipak(required.packages)

library(utils)
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE)
library(estimate)
help(estimate)

load("./Data/expression matrix/BRCA_gene_RNAseq_normalized_Tissue_Type_and_whitelists_filtered.Rdata")
RNASeq.QN.LOG2 = log(filtered.norm.RNAseqData +1, 2)
#load("~/Dropbox (TBI-Lab)/DB-LAB/Projects - Data/02 Active Projects/SDR100029_JSREP_2016_WH_Colon_Cancer_NGS/Data/Complete cohort phase/At Sidra/Data/RNASeq/Trimmed/Normalized/JSREP.clean.dataset.EDAseq.QN.HPC.Rdata")
write.table(RNASeq.QN.LOG2, sep = "\t", file = "./Data/expression matrix/BRCA_gene_RNAseq_normalized_Tissue_Type_and_whitelists_filtered_LOG2.txt", quote = FALSE)
dir.create("./Analysis/ESTIMATE", showWarnings = FALSE)

# Calculate estimate score
filterCommonGenes(input.f=paste0("./Data/expression matrix/BRCA_gene_RNAseq_normalized_Tissue_Type_and_whitelists_filtered_LOG2.txt"),
                  output.f=paste0("./Analysis/ESTIMATE/ESTIMATE.input.gct"),
                  id=c("GeneSymbol","EntrezID"))
# "Merged dataset includes 10054 genes (358 mismatched)."

estimateScore(input.ds ="./Analysis/ESTIMATE/ESTIMATE.input.gct",
              output.ds = "./Analysis/ESTIMATE/ESTIMATE.score.BRCA.TCGA.gct",
              platform= "illumina")

# [1] "1 gene set: StromalSignature  overlap= 137"
# [1] "2 gene set: ImmuneSignature  overlap= 140"

estimate.gct<-read.table("./Analysis/ESTIMATE/ESTIMATE.score.BRCA.TCGA.gct", skip = 2, header = TRUE) #skip=2 remove the first 2 lines
rownames(estimate.gct) = estimate.gct$NAME
estimate.gct$NAME = NULL
estimate.gct$Description = NULL

ESTIMATE = t(estimate.gct)
rownames(ESTIMATE) = gsub("X", "", rownames(ESTIMATE))
rownames(ESTIMATE) = gsub("\\.", "-", rownames(ESTIMATE))
save(ESTIMATE, file = "./Analysis/ESTIMATE/TCGA_BRCA_ESTIMATE_scores.Rdata")
