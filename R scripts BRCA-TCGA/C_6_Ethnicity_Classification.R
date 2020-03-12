#######
#
# Assign ethnicity to patients based on nationality
#
#######


# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA")

source("../tools/ipak.function.R")
ipak("R.utils")


# Load data
load("./Analysis/Sample_annotations.Rdata")
patient_data = read.csv("./Data/Clinical Data/patient.csv", stringsAsFactors = FALSE)

Clinical_data_ann$ethnicity = patient_data$ethnicity[match(Clinical_data_ann$bcr_patient_barcode, patient_data$bcr_patient_barcode)]

Clinical_data_ann$Ethnicity_reported = NA

Clinical_data_ann$Ethnicity_reported = paste(Clinical_data_ann$ethnicity, Clinical_data_ann$race)
Clinical_data_ann$Ethnicity_reported_simple = Clinical_data_ann$Ethnicity_reported


Clinical_data_ann[Clinical_data_ann$Ethnicity_reported == "[Not Available] [Not Available]","Ethnicity_reported_simple"] <- NA


Clinical_data_ann[Clinical_data_ann$Ethnicity_reported=="[Not Available] [Not Available]","Ethnicity_reported_simple"] <- NA
Clinical_data_ann[Clinical_data_ann$Ethnicity_reported=="[Not Evaluated] [Not Evaluated]","Ethnicity_reported_simple"] <- NA
Clinical_data_ann[Clinical_data_ann$Ethnicity_reported=="[Unknown] [Unknown]","Ethnicity_reported_simple"] <- NA
Clinical_data_ann[Clinical_data_ann$Ethnicity_reported=="[Not Evaluated] [Unknown]","Ethnicity_reported_simple"] <- NA
Clinical_data_ann[Clinical_data_ann$Ethnicity_reported=="[Unknown] [Not Evaluated]","Ethnicity_reported_simple"] <- NA


Clinical_data_ann[Clinical_data_ann$Ethnicity_reported=="[Not Available] WHITE","Ethnicity_reported_simple"] <- "WHITE"
Clinical_data_ann[Clinical_data_ann$Ethnicity_reported=="[Not Evaluated] WHITE","Ethnicity_reported_simple"] <- "WHITE"
Clinical_data_ann[Clinical_data_ann$Ethnicity_reported=="NOT HISPANIC OR LATINO WHITE","Ethnicity_reported_simple"] <- "WHITE"
Clinical_data_ann[Clinical_data_ann$Ethnicity_reported=="HISPANIC OR LATINO WHITE","Ethnicity_reported_simple"] <- "HISPANIC WHITE"
Clinical_data_ann[Clinical_data_ann$Ethnicity_reported=="[Unknown] WHITE","Ethnicity_reported_simple"] <- "WHITE"

Clinical_data_ann[Clinical_data_ann$Ethnicity_reported=="NOT HISPANIC OR LATINO BLACK OR AFRICAN AMERICAN","Ethnicity_reported_simple"] <- "BLACK"
Clinical_data_ann[Clinical_data_ann$Ethnicity_reported=="HISPANIC OR LATINO BLACK OR AFRICAN AMERICAN","Ethnicity_reported_simple"] <- "HISPANIC BLACK"
Clinical_data_ann[Clinical_data_ann$Ethnicity_reported=="[Unknown] BLACK OR AFRICAN AMERICAN","Ethnicity_reported_simple"] <- "BLACK"
Clinical_data_ann[Clinical_data_ann$Ethnicity_reported=="[Not Available] BLACK OR AFRICAN AMERICAN","Ethnicity_reported_simple"] <- "BLACK"
Clinical_data_ann[Clinical_data_ann$Ethnicity_reported=="[Not Evaluated] BLACK OR AFRICAN AMERICAN","Ethnicity_reported_simple"] <- "BLACK"

Clinical_data_ann[Clinical_data_ann$Ethnicity_reported=="NOT HISPANIC OR LATINO ASIAN","Ethnicity_reported_simple"] <- "ASIAN"
Clinical_data_ann[Clinical_data_ann$Ethnicity_reported=="HISPANIC OR LATINO ASIAN","Ethnicity_reported_simple"] <- "HISPANIC ASIAN"
Clinical_data_ann[Clinical_data_ann$Ethnicity_reported=="[Not Available] ASIAN","Ethnicity_reported_simple"] <- "ASIAN"
Clinical_data_ann[Clinical_data_ann$Ethnicity_reported=="[Not Evaluated] ASIAN","Ethnicity_reported_simple"] <- "ASIAN"
Clinical_data_ann[Clinical_data_ann$Ethnicity_reported=="[Unknown] ASIAN","Ethnicity_reported_simple"] <- "ASIAN"

Clinical_data_ann[Clinical_data_ann$Ethnicity_reported=="[Not Available] NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER","Ethnicity_reported_simple"] <- "PACIFIC ISLANDER"
Clinical_data_ann[Clinical_data_ann$Ethnicity_reported=="NOT HISPANIC OR LATINO NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER","Ethnicity_reported_simple"] <- "PACIFIC ISLANDER"
Clinical_data_ann[Clinical_data_ann$Ethnicity_reported=="HISPANIC OR LATINO NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER","Ethnicity_reported_simple"] <- "HISPANIC PACIFIC ISLANDER"

Clinical_data_ann[Clinical_data_ann$Ethnicity_reported=="[Unknown] AMERICAN INDIAN OR ALASKA NATIVE","Ethnicity_reported_simple"] <- "AMERICAN NATIVE"
Clinical_data_ann[Clinical_data_ann$Ethnicity_reported=="NOT HISPANIC OR LATINO AMERICAN INDIAN OR ALASKA NATIVE","Ethnicity_reported_simple"] <- "AMERICAN NATIVE"
Clinical_data_ann[Clinical_data_ann$Ethnicity_reported=="[Not Evaluated] AMERICAN INDIAN OR ALASKA NATIVE","Ethnicity_reported_simple"] <- "AMERICAN NATIVE"
Clinical_data_ann[Clinical_data_ann$Ethnicity_reported=="HISPANIC OR LATINO AMERICAN INDIAN OR ALASKA NATIVE","Ethnicity_reported_simple"] <- "HISPANIC AMERICAN NATIVE"

Clinical_data_ann[Clinical_data_ann$Ethnicity_reported=="HISPANIC OR LATINO [Not Evaluated]","Ethnicity_reported_simple"] <- "HISPANIC"
Clinical_data_ann[Clinical_data_ann$Ethnicity_reported=="HISPANIC OR LATINO [Unknown]","Ethnicity_reported_simple"] <- "HISPANIC"
Clinical_data_ann[Clinical_data_ann$Ethnicity_reported=="HISPANIC OR LATINO [Not Available]","Ethnicity_reported_simple"] <- "HISPANIC"

Clinical_data_ann[Clinical_data_ann$Ethnicity_reported=="NOT HISPANIC OR LATINO [Not Available]","Ethnicity_reported_simple"] <- "NOT HISPANIC"
Clinical_data_ann[Clinical_data_ann$Ethnicity_reported=="NOT HISPANIC OR LATINO [Not Evaluated]","Ethnicity_reported_simple"] <- "NOT HISPANIC"
Clinical_data_ann[Clinical_data_ann$Ethnicity_reported=="NOT HISPANIC OR LATINO [Unknown]","Ethnicity_reported_simple"] <- "NOT HISPANIC"

#further grouping of hispanic
Clinical_data_ann$Ethnicity_reported_simple <- replace(as.character(Clinical_data_ann$Ethnicity_reported_simple), Clinical_data_ann$Ethnicity_reported_simple == "HISPANIC WHITE", "HISPANIC")
Clinical_data_ann$Ethnicity_reported_simple <- replace(as.character(Clinical_data_ann$Ethnicity_reported_simple), Clinical_data_ann$Ethnicity_reported_simple == "HISPANIC ASIAN", "HISPANIC")
Clinical_data_ann$Ethnicity_reported_simple <- replace(as.character(Clinical_data_ann$Ethnicity_reported_simple), Clinical_data_ann$Ethnicity_reported_simple == "HISPANIC BLACK", "HISPANIC")
Clinical_data_ann$Ethnicity_reported_simple <- replace(as.character(Clinical_data_ann$Ethnicity_reported_simple), Clinical_data_ann$Ethnicity_reported_simple == "AMERICAN NATIVE", "OTHER")
Clinical_data_ann$Ethnicity_reported_simple <- replace(as.character(Clinical_data_ann$Ethnicity_reported_simple), Clinical_data_ann$Ethnicity_reported_simple == "NOT HISPANIC", "OTHER")

Clinical_data_ann$Ethnicity_reported_simple = tolower(Clinical_data_ann$Ethnicity_reported_simple)

Clinical_data_ann$Ethnicity_reported_simple = factor(Clinical_data_ann$Ethnicity_reported_simple, levels = c("white", "black", "asian", "hispanic", "other"))

levels(Clinical_data_ann$Ethnicity_reported_simple) = c("White", "Black", "Asian", "Hispanic", "Other")
save(Clinical_data_ann, file = "./Analysis/Sample_annotations.Rdata")



