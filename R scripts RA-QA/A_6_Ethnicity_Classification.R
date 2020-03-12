#######
#
# Assign ethnicity to patients based on nationality
#
#######


# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-RA-QA")

# Load data
load("./Analysis/Sample_annotations.Rdata")


# Ethnicity reported simple
Clinical_data_ann$Ethnicity = NA
Clinical_data_ann$Ethnicity[which(Clinical_data_ann$Nat. %in% c("Qatar", "Egyptian", "Egypt", "Iraq", "Emarati",
                                                                "Jordan", "Palestine (Ghaza Str", "Sudan"))] = "Arab"
Clinical_data_ann$Ethnicity[which(Clinical_data_ann$Nat. %in% c("India", "Philippines"))] = "Asian"
Clinical_data_ann$Ethnicity[which(Clinical_data_ann$Nat. %in% c("Irish", "Iran, Islamic Republ",
                                                                "Canada"))] = "Other"
save(Clinical_data_ann, file = "./Analysis/Sample_annotations.Rdata")