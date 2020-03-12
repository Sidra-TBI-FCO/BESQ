
#Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA")                                                                   # Setwd to location were output files have to be saved.

source("../tools/ipak.function.R")

ipak(c("ggplot2", "ggpubr", "dplyr"))

# load data
load("./Analysis/Sample_annotations.Rdata")
triple_negative = read.csv("../RNAseq-Public Data/RNASeq_subset_clinicaldata.csv", stringsAsFactors = FALSE)
patient.data = read.csv("./Data/Clinical Data/patient_data.csv", stringsAsFactors = FALSE)

# Subset to only get triple negative
Clinical_data_ann_TNBC = Clinical_data_ann[which(Clinical_data_ann$bcr_patient_barcode %in% triple_negative$X),]
Clinical_data_ann_TNBC = Clinical_data_ann_TNBC[which(Clinical_data_ann_TNBC$Assigned_Ethnicity_simplified %in% c("White", "Black")),]
Clinical_data_ann_TNBC$Assigned_Ethnicity_simplified = factor(Clinical_data_ann_TNBC$Assigned_Ethnicity_simplified, levels = c("White", "Black"))
Clinical_data_ann_TNBC = Clinical_data_ann_TNBC[-which(is.na(Clinical_data_ann_TNBC$IMS_Mathews)),]

table(Clinical_data_ann_TNBC$Assigned_Ethnicity_simplified, Clinical_data_ann_TNBC$IMS_Mathews)

DF1 <- Clinical_data_ann_TNBC %>%
  group_by(Assigned_Ethnicity_simplified, IMS_Mathews) %>%
  summarise(count=n()) %>%
  mutate(perc=count/sum(count))

colors = c("BasalHer2" = "#F4A101", "BasalLumHer2" = "#92D050",
           "BasalMyo" = "#3D64E4", "MyoLumB" = "#C720F1",
           "Lum" = "#158400","LumBasal" = "#D36733",
           "MyoLumA" = "#6F1287", "MyoLumHer2" = "#F377C6")
plot = ggplot(DF1, aes(x = Assigned_Ethnicity_simplified, y =perc*100, fill = IMS_Mathews)) + geom_bar(stat="identity") +
  labs(x = "Ethnicity", y = "Percentage", fill = "IMS_Mathews", face = "bold") +
  theme(panel.grid = element_line(linetype = "solid", colour = "white"), 
        panel.background = element_rect(fill = "white", colour = "grey"),
        axis.text = element_text(size = 19, colour = "black"),
        axis.title = element_text(size = 19, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 19, colour = "black"))+
  scale_fill_manual(values= colors)

dir.create("./Figures/C34_TNBC_Proportion_BasalMyo_in_TNBC", showWarnings = FALSE)
png("./Figures/C34_TNBC_Proportion_BasalMyo_in_TNBC/C34_TNBC_Proportion_BasalMyo_in_TNBC.png", width = 5.5, height = 5, units = "in",
    res = 600)
plot(plot)
dev.off()
