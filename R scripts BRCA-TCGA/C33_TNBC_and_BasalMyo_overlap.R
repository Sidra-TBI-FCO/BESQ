# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA")                                                                   # Setwd to location were output files have to be saved.

source("../tools/ipak.function.R")

ipak(c("ggplot2", "ggpubr", "dplyr"))

# Set parameters
Gene_set = "Selected_pathways" #"Bindea_ORIG"  #"Selected_pathways"

# Load data
load(paste0("./Analysis/ssGSEA_scores/", Gene_set, "_ES_scores.Rdata"))
load("./Analysis/Sample_annotations.Rdata")
triple_negative = read.csv("../RNAseq-Public Data/RNASeq_subset_clinicaldata.csv", stringsAsFactors = FALSE)
patient.data = read.csv("./Data/Clinical Data/patient_data.csv", stringsAsFactors = FALSE)
triple_negative = read.csv("../RNAseq-Public Data/RNASeq_subset_clinicaldata.csv", stringsAsFactors = FALSE)

Survival_data_1 = Clinical_data_ann 
Survival_data_1 = Survival_data_1[which(Survival_data_1[, "IMS"] == "Basal"),]
#Survival_data_1 = Survival_data_1[-which(Survival_data_1$Assigned_Ethnicity_simplified %in% c("Unclear", "Asian")),]
#Survival_data_1 = Survival_data_1[-which(is.na(Survival_data_1$Assigned_Ethnicity_simplified)),]
triple_neg_patients = Survival_data_1$bcr_patient_barcode[which(Survival_data_1$bcr_patient_barcode %in% triple_negative$X)]
non_TNBC_patients = Survival_data_1$bcr_patient_barcode[-which(Survival_data_1$bcr_patient_barcode %in% triple_negative$X)]

patient.data.sub = patient.data[which(patient.data$bcr_patient_barcode %in% non_TNBC_patients),]
patient.data.sub = patient.data.sub[, c("breast_carcinoma_estrogen_receptor_status", "breast_carcinoma_progesterone_receptor_status", 
                                        "lab_proc_her2_neu_immunohistochemistry_receptor_status", "lab_procedure_her2_neu_in_situ_hybrid_outcome_type")]
patient.data.sub$all = paste(patient.data.sub$breast_carcinoma_estrogen_receptor_status,
                             patient.data.sub$breast_carcinoma_progesterone_receptor_status,
                             patient.data.sub$lab_proc_her2_neu_immunohistochemistry_receptor_status,
                             patient.data.sub$lab_procedure_her2_neu_in_situ_hybrid_outcome_type)
patient.data.sub.pos = patient.data.sub[grepl("Positive", patient.data.sub$all),]

non_evaluated = nrow(patient.data.sub) - nrow(patient.data.sub.pos)
percentage_positive = nrow(patient.data.sub.pos) / (nrow(Survival_data_1) - non_evaluated)
percentage_negative = 1 - percentage_positive                                                                         
