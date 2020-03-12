
# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA")                                                                   # Setwd to location were output files have to be saved.

source("../tools/ipak.function.R")

# Load data
load("./Analysis/Sample_annotations.Rdata")
load("./Analysis/ssGSEA_scores/Selected.pathways_ES_scores.Rdata")

Clinical_data_ann$ICR_ES = ES[55,][match(Clinical_data_ann$bcr_patient_barcode,
                                                substring(colnames(ES), 1, 12))]

plot = ggplot(Clinical_data_ann, aes(x = ICRscore, y = ICR_ES)) + geom_point(size = 1) + theme_bw()

save(Clinical_data_ann, file = "./Analysis/Sample_annotations.Rdata")
