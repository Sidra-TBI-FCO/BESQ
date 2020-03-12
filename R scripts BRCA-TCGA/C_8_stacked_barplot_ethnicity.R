
#IMS Stacked bar for Each ethnicity

# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA")
source("../tools/ipak.function.R")
ipak(c("dplyr", "ggplot2"))

#Set parameters
Ethnicity_variable = "Assigned_Ethnicity_simplified" # "Ethnicity_reported_simple"
IMS_group = "IMS_Mathews" # "IMS" or "IMS_Mathews"

#load data
load("./Analysis/Sample_annotations.Rdata")

# Analysis
Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann[, Ethnicity_variable])),]
if(Ethnicity_variable == "Ethnicity_reported_simple"){
  Clinical_data_ann = Clinical_data_ann[-which(Clinical_data_ann[, Ethnicity_variable] == "Other"),]
}
Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann[, IMS_group])),]
if(Ethnicity_variable == "Assigned_Ethnicity_simplified"){
  Clinical_data_ann = Clinical_data_ann[-which(Clinical_data_ann[, Ethnicity_variable] == "Unclear"),]
}

if(Ethnicity_variable == "Ethnicity_reported_simple" & IMS_group == "IMS"){
  DF1 <- Clinical_data_ann %>%
    group_by(Ethnicity_reported_simple, IMS) %>%
    summarise(count=n()) %>%
    mutate(perc=count/sum(count))
  
  colors = c("Basal" = "#FF2500",  "Her2" = "#FFBFCB", 
             "LumA" = "#01178B","LumB" = "#ADD8E6",
             "Normal" = "#06F900")
  plot = ggplot(DF1, aes(x = Ethnicity_reported_simple, y =perc*100, fill = IMS)) + geom_bar(stat="identity") +
    labs(x = "Ethnicity", y = "Percentage", fill = "IMS", face = "bold") +
    theme(panel.grid = element_line(linetype = "solid", colour = "white"), 
          panel.background = element_rect(fill = "white", colour = "grey"),
          axis.text = element_text(size = 19, colour = "black"),
          axis.title = element_text(size = 19, colour = "black"),
          legend.title = element_blank(),
          legend.text = element_text(size = 19, colour = "black"))+
    scale_fill_manual(values= colors)
  
}
if(Ethnicity_variable == "Ethnicity_reported_simple" & IMS_group == "IMS_Mathews"){
  DF1 <- Clinical_data_ann %>%
    group_by(Ethnicity_reported_simple, IMS_Mathews) %>%
    summarise(count=n()) %>%
    mutate(perc=count/sum(count))
  
  colors = c("BasalHer2" = "#F4A101", "BasalLumHer2" = "#92D050",
             "BasalMyo" = "#3D64E4", "MyoLumB" = "#C720F1",
             "Lum" = "#158400","LumBasal" = "#D36733",
             "MyoLumA" = "#6F1287", "MyoLumHer2" = "#F377C6")
  plot = ggplot(DF1, aes(x = Ethnicity_reported_simple, y =perc*100, fill = IMS_Mathews)) + geom_bar(stat="identity") +
    labs(x = "Ethnicity", y = "Percentage", fill = "IMS_Mathews", face = "bold") +
    theme(panel.grid = element_line(linetype = "solid", colour = "white"), 
          panel.background = element_rect(fill = "white", colour = "grey"),
          axis.text = element_text(size = 19, colour = "black"),
          axis.title = element_text(size = 19, colour = "black"),
          legend.title = element_blank(),
          legend.text = element_text(size = 19, colour = "black"))+
    scale_fill_manual(values= colors)
  
}
if(Ethnicity_variable == "Assigned_Ethnicity_simplified" & IMS_group == "IMS_Mathews"){
  DF1 <- Clinical_data_ann %>%
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
  
}

dir.create("./Figures/Stacked_barplots", showWarnings = FALSE)

png(paste0("./Figures/Stacked_barplots/Stacked_barplots_v2_", IMS_group,"_by_ethnicities_", Ethnicity_variable,".png"),res=600,height=5,width=6.5,unit="in")
plot(plot)
dev.off()
