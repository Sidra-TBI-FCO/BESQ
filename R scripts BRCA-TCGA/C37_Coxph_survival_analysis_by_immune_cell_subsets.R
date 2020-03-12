
####################################################################
###
### This Script calculates 
### 
### Input data:
###
### Output data are saved as Rdata file:
#####################################################################

# Before running this script, first download TCGA assembler 2.0.3 scripts http://www.compgenome.org/TCGA-Assembler/
# Setup environment

# Script for IMS Mathews!
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA")                                                                   # Setwd to location were output files have to be saved.

source("../tools/ipak.function.R")
source("../tools/ggkm_C12.R")

required.packages = c("survival","reshape","ggplot2","plyr","Rcpp","colorspace","texreg")
required.bioconductor.packages = "survival"
ipak(required.packages)
ibiopak(required.bioconductor.packages)

# Set Parameters                                                                                            # "HML_classification" or "k3" or "k4" or "k5"
Surv_cutoff_years = 10
Outcome = "OS"
Ethnicity_variable = "Assigned_Ethnicity_simplified" #"Ethnicity_reported_simple"
Ethnicity = "All"  # "All", "White", "Black"
IMS = "All" # "BasalMyo"

# Load data
#load(paste0(code_path, "Datalists/ICR_genes.RData")) 
load("./Analysis/Sample_annotations.Rdata")
load("./Analysis/ssGSEA_scores/Selected_pathways_ES_scores.Rdata")
ES_ICR = ES

load("./Analysis/ssGSEA_scores/Bindea_ORIG_ES_scores.Rdata")

#load(paste0("./Analysis/Benefit_Cluster/", IMS,".Rdata"))
#load("./Analysis/Benefit_Cluster/significant_pathways_and_importance_survival_cox_regression_BasalMyo.Rdata")

# Create folders
dir.create("./Analysis/",showWarnings = FALSE)                                   

# Analysis
dim(Clinical_data_ann)
Survival_data_1 = Clinical_data_ann
#Survival_data_1 = Survival_data_1[-which(is.na(Survival_data_1[, Ethnicity_variable])),]
Survival_data_1 = Survival_data_1[which(Survival_data_1$bcr_patient_barcode %in% 
                                         substring(colnames(ES), 1, 12)),]

## Perform survival analysis
colnames(Survival_data_1)[which(colnames(Survival_data_1) == "death_days_to")] = "days_to_death"
colnames(Survival_data_1)[which(colnames(Survival_data_1) == "last_contact_days_to")] = "days_to_last_followup"

if(Ethnicity == "All"){
  Survival_data = Survival_data_1
}else{
  Survival_data = Survival_data_1[which(Survival_data_1[,Ethnicity_variable] == Ethnicity),]
}
if(IMS == "BasalMyo"){
  Survival_data = Survival_data[which(Survival_data$IMS_Mathews == "BasalMyo"),]
}

#TS.Surv$`TGF-beta` = (TS.Surv$`TGF-beta` - min(TS.Surv$`TGF-beta`))/(max(TS.Surv$`TGF-beta`)-min(TS.Surv$`TGF-beta`))

ES = as.data.frame(t(ES))
ES$ICR = ES_ICR["ICR genes",][match(rownames(ES), colnames(ES_ICR))] 

i=1
for (i in 1:ncol(ES)){
  col = colnames(ES)[i]
  ES[, col] = (ES[, col] - min(ES[, col]))/(max(ES[,col])-min(ES[,col]))
}

ES = as.matrix(t(ES))

results = data.frame(Cell_types = rownames(ES), p_value = 0, HR=0 , CI_lower=0, CI_upper = 0)

k=1
for (k in 1:nrow(ES)){
  cell.type = rownames(ES)[k]
  Survival_data$cell_type = ES[cell.type,][match(Survival_data$bcr_patient_barcode,
                                                 substring(colnames(ES), 1, 12))]
  Y = Surv_cutoff_years * 365
  TS.EventFree = Survival_data[Survival_data[, Outcome] == "0", c(Outcome, paste0(Outcome, ".time"), "cell_type", "ICRscore")]
  colnames(TS.EventFree) = c("Status","Time", "Group", "ICRscore")
  TS.EventFree$Time = as.numeric(as.character(TS.EventFree$Time))
  TS.EventFree$Time[TS.EventFree$Time > Y] = Y
  
  TS.EventOccured = Survival_data[Survival_data[, Outcome] == "1", c(Outcome, paste0(Outcome, ".time"), "cell_type", "ICRscore")]
  colnames(TS.EventOccured) = c("Status","Time", "Group", "ICRscore")
  TS.EventOccured$Time = as.numeric(as.character(TS.EventOccured$Time))
  TS.EventOccured$Status[which(TS.EventOccured$Time> Y)] = "0"
  TS.EventOccured$Time[TS.EventOccured$Time > Y] = Y
  
  TS.Surv = rbind (TS.EventOccured,TS.EventFree)
  TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
  TS.Surv$Status <- TS.Surv$Status == "1"
  TS.Surv = subset(TS.Surv,TS.Surv$Time > 1)                                                                                         # remove patients with less then 1 day follow up time
  
  # Cox regression
  uni_variate = coxph(formula = Surv(Time, Status) ~ Group, data = TS.Surv)
  summary = summary(uni_variate)
  mHR.extract = extract.coxph(uni_variate, include.aic = TRUE,
                              include.rsquared = TRUE, include.maxrs=TRUE,
                              include.events = TRUE, include.nobs = TRUE,
                              include.missings = TRUE, include.zph = TRUE)
  beta = coef(uni_variate)
  se   = sqrt(diag(uni_variate$var))
  p    = 1 - pchisq((beta/se)^2, 1)
  
  results$p_value[which(results$Cell_types == cell.type)] = p
  results$HR[which(results$Cell_types == cell.type)] = summary$coefficients[2]
  results$CI_lower[which(results$Cell_types == cell.type)] = summary$conf.int[3]
  results$CI_upper[which(results$Cell_types == cell.type)] = summary$conf.int[4]
 
}
dir.create("./Analysis/C37_Coxph_Bindea", showWarnings = FALSE)
save(results, file = paste0("./Analysis/C37_Coxph_Bindea/Coxph_table_", IMS, "_", Ethnicity, ".Rdata"))

