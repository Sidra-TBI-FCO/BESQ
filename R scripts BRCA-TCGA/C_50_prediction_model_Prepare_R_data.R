
# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA")                                                                   # Setwd to location were output files have to be saved.

source("../tools/ipak.function.R")

required.packages = c("ggplot2", "Hmisc")
ipak(required.packages)

# Load data
load("./Data/expression matrix/BRCA_gene_RNAseq_normalized_Tissue_Type_and_whitelists_filtered.Rdata")
RNASeq.QN.LOG2 = log(filtered.norm.RNAseqData +1, 2)
load("./Analysis/Sample_annotations.Rdata")
load("./Analysis/ssGSEA_scores/Selected_pathways_ES_scores.Rdata")

# Set parameters
group = "blacks" # "blacks" or "whites" or "alls"

# Load based on group
no_leaf = read.csv(paste0("./Analysis/Hossam_prediction_model/res_new/", group, "/noleaf_", capitalize(gsub("s", "", group)),"_trees_to_dataframe.csv"),
                   stringsAsFactors = FALSE)

#all_trees = read.csv(paste0("./Analysis/Hossam_prediction_model/res_new/", group, "/", capitalize(gsub("s", "", group)),"_trees_to_dataframe.csv"),
                     #stringsAsFactors = FALSE)


# Test subsetting of all_trees
first_12_trees = all_trees[which(all_trees$Tree %in% 0:11),]
first_12_trees_no_leaf = first_12_trees[-which(first_12_trees$Feature == "Leaf"),]
max_splits = aggregate(first_12_trees_no_leaf$Split, by = list(first_12_trees_no_leaf$Feature), FUN = max)
colnames(max_splits) = c("Pathway", "Cutoff")

# Text editing to match spelling in enrichment scores tables
max_splits$Pathway = gsub("_", " ", max_splits$Pathway)
max_splits$Pathway = gsub("\\:", "] ", max_splits$Pathway)
max_splits$Pathway = paste0("[", max_splits$Pathway)
max_splits$Pathway[which(max_splits$Pathway == "[ICRscore")] = "ICRscore"

# Save
dir.create("./Analysis/Hossam_prediction_model/cutoffs_JR", showWarnings = FALSE)

save(max_splits,file = paste0("./Analysis/Hossam_prediction_model/cutoffs_JR/", group, "_maxim_cutoff.Rdata"))
