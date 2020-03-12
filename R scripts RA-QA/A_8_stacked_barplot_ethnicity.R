
#IMS Stacked bar for Each ethnicity

# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-RA-QA")
source("../tools/ipak.function.R")
ipak("dplyr")

#load data
load("./Analysis/Sample_annotations.Rdata")

# Set parameters
IMS_group = "IMS_Mathews" # "IMS"

# Analysis
Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann$IMS)),]
#Clinical_data_ann = Clinical_data_ann[-which(Clinical_data_ann$Ethnicity == "Other"),]
Clinical_data_ann$Ethnicity[-which(Clinical_data_ann$Ethnicity == "Arab")] = "non-Arab"

if(IMS_group == "IMS"){
  DF1 <- Clinical_data_ann %>%
    group_by(Ethnicity, IMS) %>%
    summarise(count=n()) %>%
    mutate(perc=count/sum(count))
  plot = ggplot(DF1, aes(x = Ethnicity, y =perc*100, fill = IMS)) + geom_bar(stat="identity") +
    labs(x = "Ethnicity", y = "Percentage", fill = "IMS", face = "bold") +
    theme(panel.grid = element_line(linetype = "solid", colour = "white"), 
          panel.background = element_rect(fill = "white", colour = "grey"),
          axis.text = element_text(size = 19, colour = "black"),
          axis.title = element_text(size = 19, colour = "black"),
          legend.title = element_blank(),
          legend.text = element_text(size = 19, colour = "black"))+
    scale_fill_manual(values= c("Basal" = "#FF2500",  "Her2" = "#FFBFCB", 
                                "LumA" = "#01178B","LumB" = "#ADD8E6",
                                "Normal" = "#06F900"))
}

if(IMS_group == "IMS_Mathews"){
  DF1 <- Clinical_data_ann %>%
    group_by(Ethnicity, IMS_Mathews) %>%
    summarise(count=n()) %>%
    mutate(perc=count/sum(count))
  plot = ggplot(DF1, aes(x = Ethnicity, y =perc*100, fill = IMS_Mathews)) + geom_bar(stat="identity") +
    labs(x = "Ethnicity", y = "Percentage", fill = "IMS_Mathews", face = "bold") +
    theme(panel.grid = element_line(linetype = "solid", colour = "white"), 
          panel.background = element_rect(fill = "white", colour = "grey"),
          axis.text = element_text(size = 19, colour = "black"),
          axis.title = element_text(size = 19, colour = "black"),
          legend.title = element_blank(),
          legend.text = element_text(size = 19, colour = "black"))+
    scale_fill_manual(values= c("BasalHer2" = "#F4A101", "BasalLumHer2" = "#92D050",
                                "BasalMyo" = "#3D64E4", "MyoLumB" = "#C720F1",
                                "Lum" = "#158400","LumBasal" = "#D36733",
                                "MyoLumA" = "#6F1287", "MyoLumHer2" = "#F377C6"))
}
  

dir.create("./Figures/Stacked_barplots", showWarnings = FALSE)

png(paste0("./Figures/Stacked_barplots/Oct_Stacked_barplots_v2_", IMS_group ,"_by_ethnicities.png"),res=600,height=4.5,width=5,unit="in")
plot(plot)
dev.off()
