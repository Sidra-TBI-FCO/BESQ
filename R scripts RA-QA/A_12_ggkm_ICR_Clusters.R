
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
source("../tools/ggkm_A12.R")

required.packages = c("survival","reshape","ggplot2","plyr","Rcpp","colorspace","texreg")
required.bioconductor.packages = "survival"
ipak(required.packages)
ibiopak(required.bioconductor.packages)

# Set Parameters                                                                                            # "HML_classification" or "k3" or "k4" or "k5"
Surv_cutoff_years = 10
Outcome = "OS"
IMS = "All"

# Load data
#load(paste0(code_path, "Datalists/ICR_genes.RData")) 
load("./Analysis/Sample_annotations.Rdata")

# Create folders
dir.create("./Analysis/",showWarnings = FALSE)                                                                        # Create folder to save processed data (by Assembler module B)
dir.create("./Figures/Kaplan_meiers", showWarnings = FALSE)
dir.create("./Figures/Kaplan_meiers/A12_survival_ICR_by_Ethnicity", showWarnings = FALSE)

# Define parameters (based on loaded data)
Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann$HL.ICR.Cluster)),]
Clinical_data_ann$HL.ICR.Cluster = factor(Clinical_data_ann$HL.ICR.Cluster, levels = c("ICR-High","ICR-Low"))

IMS_categories = c("All", names(table(Clinical_data_ann$IMS)))
Ethnicities = c("All", "Arab", "Asian")

All_survival_analysis_data = data.frame(Ethnicity = Ethnicities, p_value = 0, HR=0 , CI_lower=0, CI_upper = 0)
N.sets = length(Ethnicities)

dim(Clinical_data_ann)
Survival_data_1 = Clinical_data_ann
Survival_data_1 = Survival_data_1[-which(Survival_data_1$Ethnicity %in% c("Other")),]
#Survival_data_1 = Survival_data_1[which(Survival_data_1$IMS == IMS),]
dim(Survival_data_1)
## Perform survival analysis
Survival_data_1$OS = Survival_data_1$Status
Survival_data_1$OS[which(Survival_data_1$OS == "Alive")] = 0
Survival_data_1$OS[which(Survival_data_1$OS == "Dead")] = 1

i = 2
for (i in 1:N.sets) {
  Ethnicity = Ethnicities[i]
  
  if(Ethnicity == "All"){
    Survival_data = Survival_data_1
  }else{
    Survival_data = Survival_data_1[which(Survival_data_1$Ethnicity == Ethnicity),]
  }
  
  Y = Surv_cutoff_years * 365
  TS.EventFree = Survival_data[Survival_data[, Outcome] == "0", c(Outcome, paste0(Outcome, ".time"), "HL.ICR.Cluster")]
  colnames(TS.EventFree) = c("Status","Time", "Group")
  TS.EventFree$Time = as.numeric(as.character(TS.EventFree$Time))
  TS.EventFree$Time[TS.EventFree$Time > Y] = Y
  
  TS.EventOccured = Survival_data[Survival_data[, Outcome] == "1", c(Outcome, paste0(Outcome, ".time"), "HL.ICR.Cluster")]
  colnames(TS.EventOccured) = c("Status","Time", "Group")
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
  
  TS.Surv[,"Group"] = factor(TS.Surv[,"Group"], levels = c("ICR-High", "ICR-Low"))
  
  # Check this!!
  ##TS.Surv[,"Group"] = relevel(TS.Surv[,"Group"], "ICR High")
  mHR = coxph(formula = msurv ~ TS.Surv[,"Group"],data = TS.Surv, subset = TS.Surv$Group %in% c("ICR-High", "ICR-Low"))
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
  x = 1
  All_survival_analysis_data[All_survival_analysis_data == Ethnicity, "p_value"] = p[x]
  All_survival_analysis_data[All_survival_analysis_data == Ethnicity, "HR"] = signif(exp(mHR.extract@coef),3)[x]
  All_survival_analysis_data[All_survival_analysis_data == Ethnicity, "CI_lower"] = CI[x,1]
  All_survival_analysis_data[All_survival_analysis_data == Ethnicity, "CI_upper"] = CI[x,2]
  
  PLOT_P = signif(p[x], digits = 3)
  PLOT_HR = round(signif(exp(mHR.extract@coef),3)[x], 3)
  PLOT_CI1 = CI[x,1]
  PLOT_CI2 = CI[x,2]
  
  # plots
  png(paste0("./Figures/Kaplan_meiers/A12_survival_ICR_by_Ethnicity/", Ethnicity, "_", IMS, "_", Outcome,
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

