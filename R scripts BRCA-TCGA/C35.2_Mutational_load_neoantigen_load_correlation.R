
#Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA")                                                                   # Setwd to location were output files have to be saved.

source("../tools/ipak.function.R")

ipak(c("ggplot2", "ggpubr", "dplyr"))

# load data
load("./Analysis/Sample_annotations.Rdata")
patient.data = read.csv("./Data/Clinical Data/patient_data.csv", stringsAsFactors = FALSE)
PanImmune_MS = read.csv("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/3_DataProcessing/External/mmc2-PanImmune_MS.csv", stringsAsFactors = FALSE)
load("~/Dropbox (TBI-Lab)/TCGA pancancer-germline/New Master Data/Master.file.new.RS.v3.Rdata")

# Set parameters
Ethnicity_variable = "Assigned_Ethnicity_simplified"
variable_of_interest = "Nonsilent.Mutation.Rate" #"SNV.Neoantigens" #"Nonsilent.Mutation.Rate"
log_transform = "log_transform"
source = "TCGA_Masterfile"

# Analysis
Clinical_data_ann = Clinical_data_ann[-which(Clinical_data_ann[, Ethnicity_variable] %in% c("Other", "Asian", "Hispanic", "Unclear")),]
Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann$HML.ICR.Cluster)),]
Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann[, Ethnicity_variable])),]
Clinical_data_ann = Clinical_data_ann[which(Clinical_data_ann$IMS_Mathews == "BasalMyo"),]

#Clinical_data_ann = Clinical_data_ann[which(Clinical_data_ann$Assigned_Ethnicity_simplified == "White"),]

if(source == "PanImmune"){
  Clinical_data_ann$Mut_load = PanImmune_MS$Nonsilent.Mutation.Rate[match(Clinical_data_ann$bcr_patient_barcode, PanImmune_MS$TCGA.Participant.Barcode)]
  Clinical_data_ann$Neoantigens = PanImmune_MS$SNV.Neoantigens[match(Clinical_data_ann$bcr_patient_barcode, PanImmune_MS$TCGA.Participant.Barcode)]
}

if(source == "TCGA_Masterfile"){
  Clinical_data_ann$Mut_load = pheno.data$Core56.DNA.Alteration_Nonsilent.Mutation.Rate[match(Clinical_data_ann$bcr_patient_barcode, pheno.data$Patient.ID)]
  Clinical_data_ann$Neoantigens = pheno.data$Core56.DNA.Alteration_SNV.Neoantigens[match(Clinical_data_ann$bcr_patient_barcode, pheno.data$Patient.ID)]
}

Clinical_data_ann$Mut_load = log10(Clinical_data_ann$Mut_load + 1)
Clinical_data_ann$Neoantigens = log10(Clinical_data_ann$Neoantigens + 1)

Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann$Mut_load)),]
Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann$Neoantigens)),]

df_plot = Clinical_data_ann

cor(df_plot$Mut_load, df_plot$Neoantigens, method = "pearson")

plot = ggplot(df_plot, aes(x = Mut_load, y = Neoantigens)) +
  geom_point(size = 0.6) +
  theme_bw() +
  xlab("Mutation Rate (log10(x+1)") +
  ylab("Neoantigen count (log10(x+1)") +
  geom_smooth(method=lm, color="blue") +
  geom_text(label = paste0("Pearson's rho = ", round(cor(df_plot$Mut_load, df_plot$Neoantigens, method = "pearson"), 2)),
            x = 0.5, y = 2.5, size = 4) +
  ggtitle("")

dev.new()
plot(plot)
