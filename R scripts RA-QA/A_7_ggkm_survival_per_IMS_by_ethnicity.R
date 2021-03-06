
####################################################################
###
### This Script calculates 
### 
### Input data:
### ("./3_DataProcessing/",download.method,"/",Cancer,"/SurvivalData/")
### Output data are saved as Rdata file:
#####################################################################

# Before running this script, first download TCGA assembler 2.0.3 scripts http://www.compgenome.org/TCGA-Assembler/
# Setup environment
# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-RA-QA")                                                                   # Setwd to location were output files have to be saved.

source("../tools/ipak.function.R")
source("../tools/ggkm_A7.R")

required.packages = c("survival","reshape","ggplot2","plyr","Rcpp","colorspace","texreg")
required.bioconductor.packages = "survival"
ipak(required.packages)
ibiopak(required.bioconductor.packages)

# Set Parameters                                                                                            # "HML_classification" or "k3" or "k4" or "k5"
Surv_cutoff_years = 10
Outcome = "OS"

# Load data
#load(paste0(code_path, "Datalists/ICR_genes.RData")) 
load("./Analysis/Sample_annotations.Rdata")

# Create folders
dir.create("./Analysis/",showWarnings = FALSE)                                                                        # Create folder to save processed data (by Assembler module B)
dir.create("./Figures/Kaplan_meiers", showWarnings = FALSE)
dir.create("./Figures/Kaplan_meiers/A7_survival_by_IMS", showWarnings = FALSE)

# Analysis
Clinical_data_ann$OS.time = as.Date(as.character(Clinical_data_ann$Last.seen.date..updated.), format="%m/%d/%Y")-
  as.Date(as.character(Clinical_data_ann$Diagnosis.Date), format="%m/%d/%Y")

# Define parameters (based on loaded data)
IMS_categories = c("All", names(table(Clinical_data_ann$IMS)))

All_survival_analysis_data = data.frame(IMS = IMS_categories, p_value = 0, HR=0 , CI_lower=0, CI_upper = 0)
N.sets = length(IMS_categories)

dim(Clinical_data_ann)
Survival_data_1 = Clinical_data_ann
#Survival_data_1 = Survival_data_1[-which(Survival_data_1$Ethnicity == "Other"),]
Survival_data_1$Ethnicity[-which(Survival_data_1$Ethnicity == "Arab")] = "non-Arab"
Survival_data_1 = Survival_data_1[-which(is.na(Survival_data_1$HL.ICR.Cluster)),]
dim(Survival_data_1)
#Survival_data_1$Ethnicity = factor(Survival_data_1$Ethnicity, levels = c("Arab", "Asian"))
Survival_data_1$Ethnicity = factor(Survival_data_1$Ethnicity, levels = c("Arab", "non-Arab"))

## Perform survival analysis
Survival_data_1$OS = Survival_data_1$Status
Survival_data_1$OS[which(Survival_data_1$OS == "Alive")] = 0
Survival_data_1$OS[which(Survival_data_1$OS == "Dead")] = 1

i = 1
for (i in 1:N.sets) {
  IMS = IMS_categories[i]
  colnames(Survival_data_1)[which(colnames(Survival_data_1) == "death_days_to")] = "days_to_death"
  colnames(Survival_data_1)[which(colnames(Survival_data_1) == "last_contact_days_to")] = "days_to_last_followup"
  
  if(IMS == "All"){
    Survival_data = Survival_data_1
  }else{
    Survival_data = Survival_data_1[which(Survival_data_1$IMS == IMS),]
  }
  
  Y = Surv_cutoff_years * 365
  TS.EventFree = Survival_data[Survival_data[, Outcome] == "0", c(Outcome, paste0(Outcome, ".time"), "Ethnicity")]
  colnames(TS.EventFree) = c("Status","Time", "Group")
  TS.EventFree$Time = as.numeric(as.character(TS.EventFree$Time))
  TS.EventFree$Time[TS.EventFree$Time > Y] = Y
  
  TS.EventOccured = Survival_data[Survival_data[, Outcome] == "1", c(Outcome, paste0(Outcome, ".time"), "Ethnicity")]
  colnames(TS.EventOccured) = c("Status","Time", "Group")
  TS.EventOccured$Time = as.numeric(as.character(TS.EventOccured$Time))
  TS.EventOccured$Status[which(TS.EventOccured$Time> Y)] = "EventFree"
  TS.EventOccured$Time[TS.EventOccured$Time > Y] = Y
  
  TS.Surv = rbind (TS.EventOccured,TS.EventFree)
  TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
  TS.Surv$Status <- TS.Surv$Status == "1"
  TS.Surv = subset(TS.Surv,TS.Surv$Time > 1)                                                                                            # remove patients with less then 1 day follow up time
  
  # survival curve
  msurv = Surv(TS.Surv$Time/30.4, TS.Surv$Status)                                                                                    # calculate the number of months
  mfit = survfit(msurv~TS.Surv$Group,conf.type = "log-log")
  
  # Calculations
  mdiff = survdiff(eval(mfit$call$formula), data = eval(mfit$call$data))
  pval = pchisq(mdiff$chisq,length(mdiff$n) - 1,lower.tail = FALSE)
  pvaltxt = ifelse(pval < 0.0001,"p < 0.0001",paste("p =", signif(pval, 3)))
  
  #TS.Surv[,"Group"] = factor(TS.Surv[,"Group"], levels = c("ICR High", "ICR Medium", "ICR Low"))
  TS.Surv[,"Group"] = factor(TS.Surv[,"Group"], levels = c("Arab", "non-Arab"))
  
  # Check this!!
  ##TS.Surv[,"Group"] = relevel(TS.Surv[,"Group"], "ICR High")
  mHR = coxph(formula = msurv ~ TS.Surv[,"Group"],data = TS.Surv, subset = TS.Surv$Group %in% c("Arab", "non-Arab"))
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
  All_survival_analysis_data[All_survival_analysis_data == IMS, "p_value"] = p[1]
  All_survival_analysis_data[All_survival_analysis_data == IMS, "HR"] = signif(exp(mHR.extract@coef),3)[1]
  All_survival_analysis_data[All_survival_analysis_data == IMS, "CI_lower"] = CI[1,1]
  All_survival_analysis_data[All_survival_analysis_data == IMS, "CI_upper"] = CI[1,2]
  
  PLOT_P = signif(p[1], digits = 3)
  PLOT_HR = round(signif(exp(mHR.extract@coef),3)[1], 3)
  PLOT_CI1 = CI[1,1]
  PLOT_CI2 = CI[1,2]
  
  # plots
  png(paste0("./Figures/Kaplan_meiers/A7_survival_by_IMS/Oct_", IMS, "_", Outcome, "_Kaplan_Meijer_by_ethnicity.png"),
      res=600, height=3,width=4,unit="in")                                                                                           # set filename
  ggkm(mfit,
       timeby=12,
       ystratalabs = levels(TS.Surv[,"Group"]),
       ystrataname = NULL,
       main= paste0(IMS, "\n"),
       xlabs = "Time in months",
       cbPalette = cbPalette,
       legend = "none")
  #PLOT_HR = PLOT_HR,
  #PLOT_P = PLOT_P,
  #PLOT_CI1 = PLOT_CI1,
  #PLOT_CI2 = PLOT_CI2)
  dev.off()
  
  #####
}

