
# Set-up environment 
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA")
source("../tools/ipak.function.R")
required.packages <- c("ComplexHeatmap", "impute", "ctc", "amap", "dendextend")
ipak(required.packages)  

# Load data
Survival_gx_boost = read.csv("./Analysis/Gxboost_trees_analysis_Raghvendra/Results/Feature_Importance_For_Survival.csv", stringsAsFactors = FALSE)
Survival_gx_boost$Features = paste("[", Survival_gx_boost$Features, sep = "")
Survival_gx_boost$Features = gsub("\\:", "] ", Survival_gx_boost$Features)
Survival_gx_boost$Features = gsub("\\_", " ", Survival_gx_boost$Features)
load("./Analysis/Sample_annotations.Rdata")
load("./Analysis/ssGSEA_scores/Selected_pathways_ES_scores.Rdata")

# Set Parameters                                                                                            # "HML_classification" or "k3" or "k4" or "k5"
Surv_cutoff_years = 10
Outcome = "OS"
Ethnicity_variable = "Assigned_Ethnicity_simplified" #"Ethnicity_reported_simple"
IMS_groups = c("BasalMyo") # c("BasalMyo", "BasalHer2")
IMS = "BasalMyo"
Combine_ICR_ML = ""

Clinical_data_ann$HML.ICR.Cluster = factor(Clinical_data_ann$HML.ICR.Cluster, levels = c("ICR High", "ICR Medium", "ICR Low"))
dim(Clinical_data_ann)
Survival_data = Clinical_data_ann
Survival_data = Survival_data[-which(is.na(Survival_data[, Ethnicity_variable])),]
#Survival_data_1 = Survival_data_1[-which(Survival_data_1[, Ethnicity_variable] %in% c("Unclear")),]
if(IMS == "All"){
}else{Survival_data = Survival_data[which(Survival_data$IMS_Mathews %in% IMS_groups),]}

Survival_data = Survival_data[which(Survival_data[,Ethnicity_variable] %in% c("White", "Black")),]
dim(Survival_data)

ES = as.data.frame(t(ES))
rownames(ES) = substring(rownames(ES), 1, 12)
ES$bcr_patient_barcode = rownames(ES)
ES = ES[Survival_data$bcr_patient_barcode,]
Survival_data = merge(Survival_data, ES, by = "bcr_patient_barcode")
colnames(Survival_data)[which(colnames(Survival_data) == Outcome)] = "Status"
colnames(Survival_data)[which(colnames(Survival_data) == paste0(Outcome, ".time"))] = "Time"
colnames(Survival_data)[which(colnames(Survival_data) == "HML.ICR.Cluster")] = "Group"

Y = Surv_cutoff_years * 365
TS.EventFree = Survival_data[Survival_data$Status == "0",]
TS.EventFree$Time = as.numeric(as.character(TS.EventFree$Time))
TS.EventFree$Time[TS.EventFree$Time > Y] = Y

TS.EventOccured = Survival_data[Survival_data$Status == "1",]
TS.EventOccured$Time = as.numeric(as.character(TS.EventOccured$Time))
TS.EventOccured$Status[which(TS.EventOccured$Time> Y)] = "0"
TS.EventOccured$Time[TS.EventOccured$Time > Y] = Y

TS.Surv = rbind (TS.EventOccured,TS.EventFree)
TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
TS.Surv$Status <- TS.Surv$Status == "1"
TS.Surv = subset(TS.Surv,TS.Surv$Time > 1)                                                                                         # remove patients with less then 1 day follow up time

# survival curve
msurv = Surv(TS.Surv$Time/30.4, TS.Surv$Status)                                                                                    # calculate the number of months
mfit = survfit(msurv~TS.Surv$Group,conf.type = "log-log")

# Calculations
mdiff = survdiff(eval(mfit$call$formula), data = eval(mfit$call$data))
pval = pchisq(mdiff$chisq,length(mdiff$n) - 1,lower.tail = FALSE)
pvaltxt = ifelse(pval < 0.0001,"p < 0.0001",paste("p =", signif(pval, 3)))

# Check this!!
##TS.Surv[,"Group"] = relevel(TS.Surv[,"Group"], "ICR High")

mHR = coxph(formula = msurv ~ TS.Surv[,"Group"],data = TS.Surv, subset = TS.Surv$Group %in% c("ICR High", "ICR Low"))
mHR.extract = extract.coxph(mHR, include.aic = TRUE,
                            include.rsquared = TRUE, include.maxrs=TRUE,
                            include.events = TRUE, include.nobs = TRUE,
                            include.missings = TRUE, include.zph = TRUE)
HRtxt = paste("Hazard-ratio =", signif(exp(mHR.extract@coef),3),"for",names(mHR$coefficients))
beta = coef(mHR)
se   = sqrt(diag(mHR$var))
p    = 1 - pchisq((beta/se)^2, 1)
CI   = confint(mHR)
CI   = round(exp(CI),2)

## extract only ICR High vs. ICR Low
if(Combine_ICR_ML == "ICR_Medium_Low_combined"){
  x = 1
}else{ x = 2 }
All_survival_analysis_data[All_survival_analysis_data == Ethnicity, "p_value"] = p[x]
All_survival_analysis_data[All_survival_analysis_data == Ethnicity, "HR"] = signif(exp(mHR.extract@coef),3)[x]
All_survival_analysis_data[All_survival_analysis_data == Ethnicity, "CI_lower"] = CI[x,1]
All_survival_analysis_data[All_survival_analysis_data == Ethnicity, "CI_upper"] = CI[x,2]

PLOT_P = signif(p[x], digits = 3)
PLOT_HR = round(signif(exp(mHR.extract@coef),3)[x], 3)
PLOT_CI1 = CI[x,1]
PLOT_CI2 = CI[x,2]

# Cox regression
pathways = colnames(TS.Surv)[45:99]
results = data.frame(pathway = pathways, HR = NA, p = NA)

i = 1
for (i in 1:length(pathways)){
  pathway = pathways[i]
  TS.Surv_inloop = TS.Surv
  colnames(TS.Surv_inloop)[which(colnames(TS.Surv_inloop) == pathway)] = "Variable"
  univariate = coxph(formula = Surv(Time, Status) ~ Variable, data = TS.Surv_inloop)
  summary = summary(univariate)
  results[which(results$pathway == pathway), 2:3] = c(summary$conf.int[1], summary$coefficients[5])
}

results$Direction = results$HR < 1
results$Direction[which(results$Direction == TRUE)] = "associated with better survival"
results$Direction[which(results$Direction == FALSE)] = "associated with worse survival"

results_sig = results[which(results$p < 0.05),]
