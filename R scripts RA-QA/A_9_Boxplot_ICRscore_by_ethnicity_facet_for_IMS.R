

# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-RA-QA")                                                                   # Setwd to location were output files have to be saved.

source("../tools/ipak.function.R")

ipak(c("ggplot2", "ggpubr"))

# Set parameters
ICR_variable = "ICR_ES" #"ICRscore" # "ICR_ES"

# Load data
load("./Analysis/Sample_annotations.Rdata")
load("./Analysis/ICR data/CTA_TNBC_BRCA_Qatar_RNASeq.tumor.Clustering.RData")

#Clinical_data_ann$ICRscore = clustering$ICR.score[match(Clinical_data_ann$Sample_ID,
#                                                        clustering$Sample_ID)]
#Clinical_data_ann$HL.ICR.Cluster = clustering$HL.ICR.Cluster[match(Clinical_data_ann$Sample_ID,
#                                                                    clustering$Sample_ID)]

#save(Clinical_data_ann, file = "./Analysis/Sample_annotations.Rdata")

Clinical_data_ann = Clinical_data_ann[-which(Clinical_data_ann$Ethnicity == "Other"),]
#Clinical_data_ann$Ethnicity[-which(Clinical_data_ann$Ethnicity == "Arab")] = "non-Arab"
Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann$ICRscore)),]

plot_df = data.frame(Patient_id = Clinical_data_ann$Sample_ID,
                     Ethnicity = Clinical_data_ann$Ethnicity,
                     IMS = Clinical_data_ann$IMS_Mathews,
                     ICR = Clinical_data_ann[,ICR_variable])

my_comparisons = list(c("Arab", "Asian"))
#my_comparisons = list(c("Arab", "non-Arab"))

plot = ggplot(plot_df, aes(y = ICR, x = Ethnicity)) +
  geom_boxplot(outlier.shape = NA, colour = "grey") +
  geom_jitter(aes(color = Ethnicity), width = 0.1, size = 0.5) +
  scale_color_manual(values = c("violet", "#7F2D80")) +
  #scale_color_manual(values = c("violet", "black")) +
  #stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") +
  facet_grid(.~IMS) +
  theme_bw() +
  ylim(c(-0.3, 0.35)) +
  theme(axis.text.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = NA),
        strip.background = element_rect(colour="black", fill=NA))


dir.create("./Figures/A9_Boxplot_ICRscore_by_ethnicity_facet_for_ethnicity", showWarnings = FALSE)
png(paste0("./Figures/A9_Boxplot_ICRscore_by_ethnicity_facet_for_ethnicity/Feb_2020_A9_Boxplot_", ICR_variable,"_per_Ethnicity_by_IMS.png"),
    width = 6, height = 3, res = 600, units = "in")
print(plot)
dev.off()
