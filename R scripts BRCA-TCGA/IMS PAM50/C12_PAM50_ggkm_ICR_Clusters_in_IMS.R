
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

# Script for IMS PAM50!
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
IMS_groups = c("Basal") # c("Basal")
IMS = "Basal" # "Basal"
Combine_ICR_ML = ""
Benefit_cluster = "non-beneficial" # "All" "beneficial" "non-beneficial"

# Load data
#load(paste0(code_path, "Datalists/ICR_genes.RData")) 
load("./Analysis/Sample_annotations.Rdata")
load(paste0("./Analysis/Benefit_Cluster/With_PAM50/", IMS,".Rdata"))

# Create folders
dir.create("./Analysis/",showWarnings = FALSE)                                                                        # Create folder to save processed data (by Assembler module B)
dir.create("./Figures/Kaplan_meiers", showWarnings = FALSE)
dir.create("./Figures/Kaplan_meiers/C12_survival_ICR_by_Ethnicity", showWarnings = FALSE)
dir.create("./Figures/Kaplan_meiers/C12_survival_ICR_by_Ethnicity/PAM50", showWarnings = FALSE)

# Analysis
Clinical_data_ann$ICR_benefit_cluster = x$cluster[match(Clinical_data_ann$bcr_patient_barcode, rownames(x))]
if(Combine_ICR_ML == "ICR_Medium_Low_combined"){
  Clinical_data_ann$HML.ICR.Cluster[which(Clinical_data_ann$HML.ICR.Cluster %in% c("ICR Medium", "ICR Low"))] = "ICR Medium-Low"
  Clinical_data_ann$HML.ICR.Cluster = factor(Clinical_data_ann$HML.ICR.Cluster, levels = c("ICR High", "ICR Medium-Low"))
}else{
  Clinical_data_ann$HML.ICR.Cluster = factor(Clinical_data_ann$HML.ICR.Cluster, levels = c("ICR High", "ICR Medium", "ICR Low"))
}

IMS_categories = c("All", names(table(Clinical_data_ann$IMS)))
Ethnicities = c("All", "White", "Black", "Asian", "Unclear")

All_survival_analysis_data = data.frame(Ethnicity = Ethnicities, p_value = 0, HR=0 , CI_lower=0, CI_upper = 0)
N.sets = length(Ethnicities)

dim(Clinical_data_ann)
Survival_data_1 = Clinical_data_ann
Survival_data_1 = Survival_data_1[-which(is.na(Survival_data_1[, Ethnicity_variable])),]
#Survival_data_1 = Survival_data_1[-which(Survival_data_1[, Ethnicity_variable] %in% c("Unclear")),]
if(IMS == "All"){
}else{Survival_data_1 = Survival_data_1[which(Survival_data_1$IMS %in% IMS_groups),]}
if(Benefit_cluster == "All"){
}else{Survival_data_1 = Survival_data_1[which(Survival_data_1$ICR_benefit_cluster == Benefit_cluster),]}
dim(Survival_data_1)
## Perform survival analysis

i = 1
for (i in 1:N.sets) {
  Ethnicity = Ethnicities[i]
  colnames(Survival_data_1)[which(colnames(Survival_data_1) == "death_days_to")] = "days_to_death"
  colnames(Survival_data_1)[which(colnames(Survival_data_1) == "last_contact_days_to")] = "days_to_last_followup"
  
  if(Ethnicity == "All"){
    Survival_data = Survival_data_1
  }else{
    Survival_data = Survival_data_1[which(Survival_data_1[,Ethnicity_variable] == Ethnicity),]
  }
  
  Y = Surv_cutoff_years * 365
  TS.EventFree = Survival_data[Survival_data[, Outcome] == "0", c(Outcome, paste0(Outcome, ".time"), "HML.ICR.Cluster", "ICRscore")]
  colnames(TS.EventFree) = c("Status","Time", "Group", "ICRscore")
  TS.EventFree$Time = as.numeric(as.character(TS.EventFree$Time))
  TS.EventFree$Time[TS.EventFree$Time > Y] = Y
  
  TS.EventOccured = Survival_data[Survival_data[, Outcome] == "1", c(Outcome, paste0(Outcome, ".time"), "HML.ICR.Cluster", "ICRscore")]
  colnames(TS.EventOccured) = c("Status","Time", "Group", "ICRscore")
  TS.EventOccured$Time = as.numeric(as.character(TS.EventOccured$Time))
  TS.EventOccured$Status[which(TS.EventOccured$Time> Y)] = "EventFree"
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
  
  if(Combine_ICR_ML == "ICR_Medium_Low_combined"){
    TS.Surv[,"Group"] = factor(TS.Surv[,"Group"], levels = c("ICR High", "ICR Medium-Low"))
  }else{
    TS.Surv[,"Group"] = factor(TS.Surv[,"Group"], levels = c("ICR High", "ICR Medium", "ICR Low"))
  }
  
  # Check this!!
  ##TS.Surv[,"Group"] = relevel(TS.Surv[,"Group"], "ICR High")
  if(Combine_ICR_ML == "ICR_Medium_Low_combined"){
    mHR = coxph(formula = msurv ~ TS.Surv[,"Group"],data = TS.Surv, subset = TS.Surv$Group %in% c("ICR High", "ICR Medium-Low"))
  }else{
    mHR = coxph(formula = msurv ~ TS.Surv[,"Group"],data = TS.Surv, subset = TS.Surv$Group %in% c("ICR High", "ICR Low"))
  }
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
  uni_variate_ICRscore = coxph(formula = Surv(Time, Status) ~ ICRscore, data = TS.Surv)
  summary(uni_variate_ICRscore)
  
  # plots
  png(paste0("./Figures/Kaplan_meiers/C12_survival_ICR_by_Ethnicity/PAM50/", Ethnicity, "_", IMS, "_", Outcome, "_",
             Combine_ICR_ML, "_In_ICR_benefit_cluster_", Benefit_cluster, 
             "_Kaplan_Meijer_by_ICR.png"),
      res=600, height=3,width=4,unit="in")                                                                                           # set filename
  ggkm(mfit,
       timeby=12,
       ystratalabs = levels(TS.Surv[,"Group"]),
       ystrataname = NULL,
       main= paste0(Ethnicity, "\n"),
       xlabs = "Time in months",
       cbPalette = cbPalette)
  #PLOT_HR = PLOT_HR,
  #PLOT_P = PLOT_P,
  #PLOT_CI1 = PLOT_CI1,
  #PLOT_CI2 = PLOT_CI2)
  dev.off()
  
  #####
}

