#################################################################
###
### 
#################################################################

# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-RA-QA") 

source("~/Dropbox (TBI-Lab)/DB-LAB/bioinformatics tools/R scripts/ipak.function.R")

required.packages <- c("corrplot", "stringr")
ipak(required.packages)

library(utils)
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE)
library(estimate)
#help(estimate)

load("./Data/expression matrix/CTA_TNBC_BRCA_Qatar_RNASeq.tumor.tissue.filtered.EDAseq.QN.Rdata")
write.table(RNASeq.QN.LOG2, sep = "\t", file = "./Data/expression matrix/RNASeq_QN_LOG2.txt", quote = FALSE)
dir.create("./Analysis/ESTIMATE", showWarnings = FALSE)

# Calculate estimate score
filterCommonGenes(input.f=paste0("./Data/expression matrix/RNASeq_QN_LOG2.txt"),
                  output.f=paste0("./Analysis/ESTIMATE/ESTIMATE.input.gct"),
                  id=c("GeneSymbol","EntrezID"))
# "Merged dataset includes 10054 genes (358 mismatched)."

estimateScore(input.ds ="./Analysis/ESTIMATE/ESTIMATE.input.gct",
              output.ds = "./Analysis/ESTIMATE/ESTIMATE.score.RA.QA.gct",
              platform= "illumina")

# [1] "1 gene set: StromalSignature  overlap= 138"
# [1] "2 gene set: ImmuneSignature  overlap= 141"

estimate.gct<-read.table("./Analysis/ESTIMATE/ESTIMATE.score.RA.QA.gct", skip = 2, header = TRUE) #skip=2 remove the first 2 lines
rownames(estimate.gct) = estimate.gct$NAME
estimate.gct$NAME = NULL
estimate.gct$Description = NULL

ESTIMATE = t(estimate.gct)
save(ESTIMATE, file = "./Analysis/ESTIMATE/TCGA_BRCA_ESTIMATE_scores.Rdata")
