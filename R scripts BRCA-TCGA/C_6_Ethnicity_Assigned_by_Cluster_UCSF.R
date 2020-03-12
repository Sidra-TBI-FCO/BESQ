#######
#
# Assigned ethnicity based on DNA profile data
#
#######

# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA")

source("../tools/ipak.function.R")
ipak("dplyr")

# Load data
load("./Analysis/Sample_annotations.Rdata")
UCSF_df = read.csv("./Data/UCSF_Ancestry_Calls.csv", stringsAsFactors = FALSE)

Clinical_data_ann$Assigned_Ethicity = UCSF_df$pam.ancestry.cluster[match(Clinical_data_ann$bcr_patient_barcode,
                                                                         UCSF_df$Patient_ID)]

DF1 = Clinical_data_ann %>%
  group_by(Ethnicity_reported_simple, Assigned_Ethicity) %>%
  summarise(count=n()) %>%
  mutate(perc=count/sum(count))

Clinical_data_ann$Assigned_Ethicity_simplified = "Unclear"

Clinical_data_ann$Assigned_Ethicity_simplified[which(Clinical_data_ann$Ethnicity_reported_simple == "Asian" &
                                                       is.na(Clinical_data_ann$Assigned_Ethicity))] = "Asian"
Clinical_data_ann$Assigned_Ethicity_simplified[which(Clinical_data_ann$Ethnicity_reported_simple == "Asian" &
                                                      Clinical_data_ann$Assigned_Ethicity == "ASIAN")] = "Asian"

Clinical_data_ann$Assigned_Ethicity_simplified[which(Clinical_data_ann$Ethnicity_reported_simple == "Black" &
                                                       is.na(Clinical_data_ann$Assigned_Ethicity))] = "Black"
Clinical_data_ann$Assigned_Ethicity_simplified[which(Clinical_data_ann$Ethnicity_reported_simple == "Black" &
                                                       Clinical_data_ann$Assigned_Ethicity == "AFR")] = "Black"
Clinical_data_ann$Assigned_Ethicity_simplified[which(is.na(Clinical_data_ann$Ethnicity_reported_simple) &
                                                       Clinical_data_ann$Assigned_Ethicity == "AFR")] = "Black"

Clinical_data_ann$Assigned_Ethicity_simplified[which(is.na(Clinical_data_ann$Ethnicity_reported_simple) &
                                                       is.na(Clinical_data_ann$Assigned_Ethicity))] = NA

Clinical_data_ann$Assigned_Ethicity_simplified[which(Clinical_data_ann$Ethnicity_reported_simple == "Hispanic" &
                                                       Clinical_data_ann$Assigned_Ethicity == "AFR")] = "Black"
Clinical_data_ann$Assigned_Ethicity_simplified[which(Clinical_data_ann$Ethnicity_reported_simple == "Hispanic" &
                                                       Clinical_data_ann$Assigned_Ethicity == "ASIAN")] = "Asian"
Clinical_data_ann$Assigned_Ethicity_simplified[which(Clinical_data_ann$Ethnicity_reported_simple == "Hispanic" &
                                                       Clinical_data_ann$Assigned_Ethicity == "EUR")] = "White"

Clinical_data_ann$Assigned_Ethicity_simplified[which(Clinical_data_ann$Ethnicity_reported_simple == "White" &
                                                       Clinical_data_ann$Assigned_Ethicity == "EUR")] = "White"
Clinical_data_ann$Assigned_Ethicity_simplified[which(is.na(Clinical_data_ann$Ethnicity_reported_simple) &
                                                       Clinical_data_ann$Assigned_Ethicity == "EUR")] = "White"
Clinical_data_ann$Assigned_Ethicity_simplified[which(Clinical_data_ann$Ethnicity_reported_simple == "White" &
                                                       is.na(Clinical_data_ann$Assigned_Ethicity))] = "White"

table(Clinical_data_ann$Assigned_Ethicity_simplified, exclude = NULL)

# Correct typo
colnames(Clinical_data_ann)[43] = "Assigned_Ethnicity"
colnames(Clinical_data_ann)[44] = "Assigned_Ethnicity_simplified"

Clinical_data_ann$Assigned_Ethnicity_simplified = factor(Clinical_data_ann$Assigned_Ethnicity_simplified,
                                                         levels = c("White", "Black", "Asian", "Unclear"))

save(Clinical_data_ann, file = "./Analysis/Sample_annotations.Rdata")

