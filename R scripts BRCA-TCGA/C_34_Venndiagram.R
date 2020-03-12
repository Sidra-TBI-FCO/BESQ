

# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA")                                                                   # Setwd to location were output files have to be saved.

source("../tools/ipak.function.R")

required.packages = c("VennDiagram")
ipak(required.packages)

# Load data
load("../tools/Selected.pathways.3.4.RData")
load("./Analysis/t_test_between_ethnicities/C23_significant_pathways_BasalMyo_All.Rdata")
black = read.csv("./Analysis/Hossam_prediction_model/Final_SHAP_p_value_and_cutoffs/black.csv", stringsAsFactors = FALSE)
white = read.csv("./Analysis/Hossam_prediction_model/Final_SHAP_p_value_and_cutoffs/white.csv", stringsAsFactors = FALSE)
black_bindea = read.csv("./Analysis/Hossam_prediction_model/Final_SHAP_p_value_and_cutoffs/black_bindea.csv", stringsAsFactors = FALSE)
white_bindea = read.csv("./Analysis/Hossam_prediction_model/Final_SHAP_p_value_and_cutoffs/white_bindea.csv", stringsAsFactors = FALSE)

#SHAP_Black = black$Feature[which(black$Pval < 0.06)]
SHAP_Black = black$Feature[1:10]
SHAP_Black = c("[HM] PI3K Akt mTOR signaling", 
               "[LM] Proliferation",
               "[HM] G2M checkpoint",
               "[IPA] PI3K AKT Signaling",
               "[IPA] Telomere Extension by Telomerase",
               "[IPA] AMPK Signaling",
               "[IPA] ERK5 Signaling",
               "[IPA] ErbB Signaling",
               "[TBI] Barrier genes",
               "[HM] UV response down")
#SHAP_White = white$Features[which(white$Pval < 0.06)]
SHAP_White = white$Features[1:10]
SHAP_White = c("[TBI] Barrier genes", 
               "[HM] Reactive oxigen species pathway",
               "[IPA] EGF Signaling",
               "[HM] DNA repair",
               "[HM] Hedgehog signaling",
               "[IPA] UVC Induced MAPK Signaling",
               "[IPA] VEGF Signaling",
               "[IPA] AMPK Signaling",
               "[IPA] Estrogen Dependent Breast Cancer Signaling",
               "[HM] UV response up")



#venn.plot <- venn.diagram(list(significant_pathways,SHAP_All, SHAP_Black, SHAP_White),NULL,fill=c("orange", "white", "white", "white"), alpha=c(0.3,0.3, 0.3, 0.3), 
 #                         cex = 2, cat.fontface=8, category.names=c("Different", "SHAP All", "SHAP Black", "SHAP White"))

venn.plot <- venn.diagram(list(significant_pathways, SHAP_Black, SHAP_White),NULL,fill=c("orange", "white", "white"), alpha=c(0.3,0.3, 0.3), 
                         cex = 2, cat.fontface=8, font.family = "arial", category.names=c("","", ""), scaled = TRUE)

#venn.plot <- venn.diagram(list(significant_pathways, SHAP_Black_pos, SHAP_Black_neg, SHAP_White_pos, SHAP_White_neg),NULL,fill=c("orange", "white", "white", "white", "white"), alpha=c(0.3,0.3, 0.3, 0.3, 0.3), 
 #                         cex = 2, cat.fontface=8, category.names=c("Different","SHAP Black positive", "SHAP Black negative",
 #                                                                   "SHAP White positive", "SHAP White negative"))

intersect(intersect(significant_pathways, SHAP_All),intersect(SHAP_Black, SHAP_White))
intersect(intersect(significant_pathways, SHAP_White), SHAP_Black)

intersect(SHAP_White, significant_pathways)

intersect(significant_pathways, SHAP_All)
dev.new()
grid.draw(venn.plot)


## Old versions

SHAP_All = c("[TBI] Barrier genes", "[IPA] EGF Signaling", "[IPA] VEGF Signaling",
             "[HM] G2M checkpoint", "[IPA] ERK5 Signaling", "[LM] Proliferation",
             "[HM] DNA repair", "[HM] TGF beta signaling", "[TPW] Hypoxia/Adenosine Immune Cell Suppression",
             "[TPW] Immunogenic Cell Death (ICD)", "[HM] PI3K Akt mTOR signaling", "[HM] Hypoxia")

SHAP_White = c("[IPA] UVC Induced MAPK Signaling", "[TBI] Barrier genes", "[IPA] HER 2 Signaling in Breast Cancer",
               "[HM] Reactive oxigen species pathway", "[IPA] VEGF Signaling", "[IPA] Estrogen Dependent Breast Cancer Signaling",
               "[IPA] AMPK Signaling", "[HM] PI3K Akt mTOR signaling", "[TBI] Phopholipase", "[IPA] Myc Mediated Apoptosis Signaling",
               "[HM] mTORC1 signaling", "[IPA] EGF Signaling", "[HM] Oxidative phosphorylation", "[HM] UV response down", 
               "[TPW] Immunogenic Cell Death (ICD)", "[LM] Proliferation")

SHAP_Black = c("[IPA] Telomere Extension by Telomerase", "[HM] PI3K Akt mTOR signaling", 
               "[LM] Proliferation", "[HM] Wnt beta catenin signaling", "[TBI] Barrier genes",
               "[IPA] AMPK Signaling",
               "[IPA] PI3K AKT Signaling", "[HM] Angiogenesis", "[IPA] ErbB Signaling")

SHAP_White_pos = c("[IPA] HER 2 Signaling in Breast Cancer", "[IPA] VEGF Signaling",
                   "[HM] PI3K Akt mTOR signaling", "[TBI] Phopholipase", "[IPA] Myc Mediated Apoptosis Signaling",
                   "[HM] mTORC1 signaling", "[HM] UV response down", "[LM] Proliferation")
SHAP_White_neg = c("[IPA] UVC Induced MAPK Signaling", "[TBI] Barrier genes", "[HM] Reactive oxigen species pathway",
                   "[IPA] Estrogen Dependent Breast Cancer Signaling", "[IPA] AMPK Signaling", "[IPA] EGF Signaling",
                   "[HM] Oxidative phosphorylation", "[TPW] Immunogenic Cell Death (ICD)")

SHAP_Black_pos = c("[HM] PI3K Akt mTOR signaling", "[LM] Proliferation", "[IPA] AMPK Signaling",
                   "[IPA] PI3K AKT Signaling", "[HM] Angiogenesis")

SHAP_Black_neg = c("[IPA] Telomere Extension by Telomerase", "[HM] Wnt beta catenin signaling",
                   "[TBI] Barrier genes", "[IPA] ErbB Signaling")

