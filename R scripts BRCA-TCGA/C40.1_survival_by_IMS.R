
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
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA")                                                                   # Setwd to location were output files have to be saved.

source("../tools/ipak.function.R")
source("~/Dropbox (TBI-Lab)/DB-LAB/bioinformatics tools/R scripts/ggkm_Jessica_Pancancer.R")

required.packages = c("survival","reshape","ggplot2","plyr","Rcpp","colorspace","texreg")
required.bioconductor.packages = "survival"
ipak(required.packages)
ibiopak(required.bioconductor.packages)

# Set Parameters                                                                                            # "HML_classification" or "k3" or "k4" or "k5"
Surv_cutoff_years = 10
Outcome = "OS"
Ethnicity_variable = "Assigned_Ethnicity_simplified" #"Assigned_Ethnicity_simplified" #"Ethnicity_reported_simple"
IMS_group = "IMS_Mathews" #"IMS_Mathews" or "IMS"

# Load data
#load(paste0(code_path, "Datalists/ICR_genes.RData")) 
load("./Analysis/Sample_annotations.Rdata")

# Create folders
dir.create("./Analysis/",showWarnings = FALSE)                                                                        # Create folder to save processed data (by Assembler module B)
dir.create("./Figures/Kaplan_meiers", showWarnings = FALSE)
dir.create("./Figures/Kaplan_meiers/C40.1_survival_by_IMS", showWarnings = FALSE)

# Define parameters (based on loaded data)
IMS_categories = names(table(Clinical_data_ann[, IMS_group]))

dim(Clinical_data_ann)
Survival_data_1 = Clinical_data_ann
Survival_data_1 = Survival_data_1[-which(is.na(Survival_data_1[, Ethnicity_variable])),]
if(Ethnicity_variable == "Ethnicity_reported_simple"){
  Survival_data_1 = Survival_data_1[-which(Survival_data_1[, Ethnicity_variable] %in% c("Other")),]
  Survival_data_1[, Ethnicity_variable] = factor(Survival_data_1[, Ethnicity_variable], levels = c("White", "Black", "Asian", "Hispanic"))
}
if(Ethnicity_variable == "Assigned_Ethnicity_simplified"){
  Survival_data_1[, Ethnicity_variable] = factor(Survival_data_1[, Ethnicity_variable], levels = c("White", "Black", "Asian", "Unclear"))
}
dim(Survival_data_1)

if(IMS_group == "IMS_Mathews"){
  Survival_data_1 = Survival_data_1[-which(is.na(Survival_data_1$IMS_Mathews)),]
}

## Perform survival analysis

Survival_data = Survival_data_1

Y = Surv_cutoff_years * 365
TS.EventFree = Survival_data[Survival_data[, Outcome] == "0", c(Outcome, paste0(Outcome, ".time"), IMS_group)]
colnames(TS.EventFree) = c("Status","Time", "Group")
TS.EventFree$Time = as.numeric(as.character(TS.EventFree$Time))
TS.EventFree$Time[TS.EventFree$Time > Y] = Y

TS.EventOccured = Survival_data[Survival_data[, Outcome] == "1", c(Outcome, paste0(Outcome, ".time"), IMS_group)]
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

#TS.Surv[,"Group"] = factor(TS.Surv[,"Group"], levels = c("ICR High", "ICR Medium", "ICR Low"))

# Check this!!
##TS.Surv[,"Group"] = relevel(TS.Surv[,"Group"], "ICR High")
mHR = coxph(formula = msurv ~ TS.Surv[,"Group"],data = TS.Surv, subset = TS.Surv$Group %in% c("BasalMyo", "Lum"))
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

#PLOT_P = signif(p[1], digits = 3)
#PLOT_HR = round(signif(exp(mHR.extract@coef),3)[1], 3)
#PLOT_CI1 = CI[1,1]
#PLOT_CI2 = CI[1,2]

# plots
png(paste0("./Figures/Kaplan_meiers/C40.1_survival_by_IMS/", Ethnicity_variable,"_",IMS_group, "_", Outcome, "_Kaplan_Meijer_with_at_risk_table.png"),
    res=600, height=11,width=12,unit="in")                                                                                           # set filename
ggkm(mfit,
     timeby=12,
     ystratalabs = levels(TS.Surv[,"Group"]),
     ystrataname = NULL,
     main= paste0("Survival curve across TDA subtypes in breast cancer"),
     xlabs = "Time in months",
     palette = c("#F4A101", "#92D050", "#3D64E4", "#C720F1", 
                 "#158400","#D36733", "#6F1287", "#F377C6")
     )
dev.off()

