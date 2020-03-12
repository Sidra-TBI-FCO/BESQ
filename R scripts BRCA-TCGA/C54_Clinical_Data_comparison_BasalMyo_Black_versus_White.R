
## Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA") 
source("../tools/ipak.function.R")
ipak(c("dplyr", "ggplot2"))

# Load data
load("./Analysis/Sample_annotations.Rdata")

#Set parameters
IMS = "All"

if(IMS == "All"){Clinical_data_ann_BM = Clinical_data_ann}else{
  Clinical_data_ann_BM = Clinical_data_ann[which(Clinical_data_ann$IMS_Mathews == IMS),]
}

BM_BW = Clinical_data_ann_BM[which(Clinical_data_ann_BM$Assigned_Ethnicity_simplified %in% c("White", "Black")),]

BM_BW$ajcc_pathologic_tumor_stage[which(BM_BW$ajcc_pathologic_tumor_stage %in% c("Stage I", "Stage IA", "Stage IB"))] = "Stage I"
BM_BW$ajcc_pathologic_tumor_stage[which(BM_BW$ajcc_pathologic_tumor_stage %in% c("Stage II", "Stage IIA", "Stage IIB"))] = "Stage II"
BM_BW$ajcc_pathologic_tumor_stage[which(BM_BW$ajcc_pathologic_tumor_stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC"))] = "Stage III"
BM_BW$ajcc_pathologic_tumor_stage[which(BM_BW$ajcc_pathologic_tumor_stage %in% c("Stage IV"))] = "Stage IV"
BM_BW$ajcc_pathologic_tumor_stage[which(BM_BW$ajcc_pathologic_tumor_stage %in% c("[Discrepancy]", "[Not Available]", "Stage X"))] = NA

table(BM_BW$ajcc_pathologic_tumor_stage, BM_BW$Assigned_Ethnicity_simplified)

colors = c("Stage I" = "green", "Stage II" = "blue",
           "Stage III" = "orange", "Stage IV" = "red")

BM_BW = BM_BW[-which(is.na(BM_BW$ajcc_pathologic_tumor_stage)),]

DF1 <- BM_BW %>%
  group_by(Assigned_Ethnicity_simplified, ajcc_pathologic_tumor_stage) %>%
  summarise(count=n()) %>%
  mutate(perc=count/sum(count))

plot = ggplot(DF1, aes(x = Assigned_Ethnicity_simplified, y =perc*100, fill = ajcc_pathologic_tumor_stage)) + geom_bar(stat="identity") +
  labs(x = "Ethnicity", y = "Percentage", fill = "ajcc_pathologic_tumor_stage", face = "bold") +
  theme(panel.grid = element_line(linetype = "solid", colour = "white"), 
        panel.background = element_rect(fill = "white", colour = "grey"),
        axis.text = element_text(size = 19, colour = "black"),
        axis.title = element_text(size = 19, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 19, colour = "black"))+
  scale_fill_manual(values= colors)

dev.new()
plot(plot)


BM_BW$Assigned_Ethnicity_simplified = as.character(BM_BW$Assigned_Ethnicity_simplified)
BM_BW$Stage[which(BM_BW$ajcc_pathologic_tumor_stage %in% c("Stage I", "Stage II"))] = "Stage I or II"
BM_BW$Stage[which(BM_BW$ajcc_pathologic_tumor_stage %in% c("Stage III", "Stage IV"))] = "Stage III or IV"
tbl = table(BM_BW$Assigned_Ethnicity_simplified, BM_BW$Stage)

chisq.test(tbl)
