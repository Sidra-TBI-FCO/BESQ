
## Set environment
rm(list = ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA")
source("../tools/ipak.function.R")
ipak(c("dplyr", "ggplot2"))

# Load data
Mathews_PAM50 = read.csv("./Analysis/PAM50/Mathews_classification.csv", stringsAsFactors = FALSE)
load("./Analysis/Sample_annotations.Rdata")

# Add PAM50 Mathews to Clinical data
Clinical_data_ann$IMS_Mathews = Mathews_PAM50$classes[match(Clinical_data_ann$bcr_patient_barcode,
                                                            substring(Mathews_PAM50$X, 1, 12))]

Clinical_data_ann$IMS_Mathews = factor(Clinical_data_ann$IMS_Mathews,
                                       levels = c("BasalHer2", "BasalMyo", "BasalLumHer2",
                                                  "Lum", "LumBasal", "MyoLumA", "MyoLumB",
                                                  "MyoLumHer2"))

test = group_by(Clinical_data_ann, IMS_Mathews, IMS) %>%
  summarize(count = n())

save(Clinical_data_ann ,file ="./Analysis/Sample_annotations.Rdata")

Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann$IMS)),]

DF1 <- Clinical_data_ann %>%
  group_by(IMS, IMS_Mathews) %>%
  summarise(count=n()) %>%
  mutate(perc=count/sum(count))

png("./Figures/Stacked_barplots/Stacked_barplots_IMS_Mathews_by_IMS.png",res=600,height=4,width=5,unit="in")
plot = ggplot(DF1, aes(x = IMS, y =perc*100, fill = IMS_Mathews)) + geom_bar(stat="identity") +
  labs(x = "IMS", y = "Percentage", fill = "IMS_Mathews", face = "bold") +
  theme(panel.grid = element_line(linetype = "solid", colour = "white"), 
        panel.background = element_rect(fill = "white", colour = "grey"),
        axis.text = element_text(size = 19, colour = "black", angle = 90),
        axis.title = element_text(size = 19, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 19, colour = "black"))+
  scale_fill_manual(values= c("BasalHer2" = "#F4A101",  "BasalLumHer2" = "#92D050",
                              "BasalMyo" = "#3D64E4", "MyoLumB" = "#C720F1",
                              "Lum" = "#158400","LumBasal" = "#D36733",
                              "MyoLumA" = "#6F1287", "MyoLumHer2" = "#F377C6"))

plot(plot)
dev.off()
