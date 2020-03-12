
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
load("./Data/expression matrix/BRCA_gene_RNAseq_normalized_Tissue_Type_and_whitelists_filtered.Rdata")

# Set parameters
Ethnicity_variable = "Assigned_Ethnicity_simplified"
variable_of_interest = "Mut_load" #"SNV.Neoantigens" #"Mut_load"
log_transform = ""
source = "TCGA_Masterfile"

# Analysis
Clinical_data_ann = Clinical_data_ann[which(Clinical_data_ann$bcr_patient_barcode %in% substring(colnames(filtered.norm.RNAseqData), 1, 12)),]

Clinical_data_ann = Clinical_data_ann[-which(Clinical_data_ann[, Ethnicity_variable] %in% c("Other", "Asian", "Hispanic", "Unclear")),]
Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann$HML.ICR.Cluster)),]
Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann[, Ethnicity_variable])),]
Clinical_data_ann = Clinical_data_ann[which(Clinical_data_ann$IMS_Mathews == "BasalMyo"),]

if(source == "PanImmune"){
  Clinical_data_ann$Mut_load = PanImmune_MS$Nonsilent.Mutation.Rate[match(Clinical_data_ann$bcr_patient_barcode, PanImmune_MS$TCGA.Participant.Barcode)]
  Clinical_data_ann$Neoantigens = PanImmune_MS$SNV.Neoantigens[match(Clinical_data_ann$bcr_patient_barcode, PanImmune_MS$TCGA.Participant.Barcode)]
  Clinical_data_ann$Aneuploidy_score = PanImmune_MS$Aneuploidy.Score[match(Clinical_data_ann$bcr_patient_barcode, PanImmune_MS$TCGA.Participant.Barcode)]
}

if(source == "TCGA_Masterfile"){
  Clinical_data_ann$Mut_load = pheno.data$Core56.DNA.Alteration_Nonsilent.Mutation.Rate[match(Clinical_data_ann$bcr_patient_barcode, pheno.data$Patient.ID)]
  Clinical_data_ann$Neoantigens = pheno.data$Core56.DNA.Alteration_SNV.Neoantigens[match(Clinical_data_ann$bcr_patient_barcode, pheno.data$Patient.ID)]
  Clinical_data_ann$Wolf_MHC1.21978456 = pheno.data$Sigs160.Wolf_MHC1.21978456[match(Clinical_data_ann$bcr_patient_barcode, pheno.data$Patient.ID)]
}

Clinical_data_ann$Mut_load = log10(Clinical_data_ann$Mut_load + 1)
Clinical_data_ann$Neoantigens = log10(Clinical_data_ann$Neoantigens + 1)

#Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann$Mut_load)),]
#Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann$Neoantigens)),]
#Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann$Wolf_MHC1.21978456)),]

df_plot = data.frame(sample_barcodes = Clinical_data_ann$bcr_patient_barcode,
                     Ethnicity = Clinical_data_ann$Assigned_Ethnicity_simplified,
                     variable = Clinical_data_ann[, variable_of_interest],
                     stringsAsFactors = FALSE)

df_plot$variable = as.numeric(df_plot$variable)

## For log transformation of y-axis if not transformed change 0 to 0.0001
df_plot$variable[which(df_plot$variable == 0)] = 0.0001
#df_plot = df_plot[-which(df_plot$variable == 0),]

#df_plot = df_plot[-which(df_plot$variable > 1),]
if(variable_of_interest == "Mut_load"){
  ylab = paste0("Nonsilent mutation rate \n (log10(x+1))")
}
if(variable_of_interest == "Neoantigens"){
  ylab = paste0("SNV Neoantigens \n (log10(x+1))")
}
if(variable_of_interest == "Aneuploidy_score"){
  ylab = "Aneuploidy score"
}
if(variable_of_interest == "Wolf_MHC1.21978456"){
  ylab = "Wolf_MHC1.21978456"
}

plot = ggplot(df_plot, aes(x=Ethnicity, y = variable)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 0.3) +
  ylab(ylab) +
  theme_bw() +
  stat_compare_means(method = "t.test", comparisons = list(c("White", "Black"))) +
  theme(axis.text.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.title.x = element_text(colour = "black", size = 13),
        axis.title.y = element_text(colour = "black", size = 13))

dir.create("./Figures/C35_Mutational_load", showWarnings = FALSE)
png(paste0("./Figures/C35_Mutational_load/v2_", source, "_", variable_of_interest, "_in_BasalMyo_by_ethnicity.png"),
    width = 3, height = 4, units = "in", res = 600)
plot(plot)
dev.off()
