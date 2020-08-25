
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
source("../tools/ggkm_v4.R")

required.packages = c("survival","reshape","ggplot2","plyr","Rcpp","colorspace","texreg", "survminer")
required.bioconductor.packages = "survival"
ipak(required.packages)
ibiopak(required.bioconductor.packages)

# Set Parameters                                                                                            # "HML_classification" or "k3" or "k4" or "k5"
Surv_cutoff_years = 5
Outcome = "DSS" # "DSS" "OS" "DFI" "PFI"
Ethnicity_variable = "Assigned_Ethnicity_simplified" #"Ethnicity_reported_simple"
cbPalette = c("#066C3C", "#EF8C12")
ICR_cluster = "All" # "All", "ICR High", "ICR Medium", "ICR Low", "ICR Medium-Low"
IMS_group = "IMS" # "IMS", "IMS_Mathews"
Stages = ""
TNBC = ""
Stage_translation = "ordinal"

# Load data
#load(paste0(code_path, "Datalists/ICR_genes.RData")) 
load("./Analysis/Sample_annotations.Rdata")
triple_negative = read.csv("../RNAseq-Public Data/RNASeq_subset_clinicaldata.csv", stringsAsFactors = FALSE)
load("./Figures/Kaplan_meiers/C39_survival_by_neoantigen_mutational_load_ratio/NeoAgs_Observed_Expected_patients.Rdata")
load("./Figures/Kaplan_meiers/C39.5_Aneuploidy_score/aneuploidy_patients.Rdata")

# Create folders
dir.create("./Analysis/",showWarnings = FALSE)                                                                        # Create folder to save processed data (by Assembler module B)
dir.create("./Figures/Kaplan_meiers", showWarnings = FALSE)
dir.create("./Figures/Kaplan_meiers/C7_survival_by_IMS", showWarnings = FALSE)

# Define parameters (based on loaded data)
IMS_categories = c("All", names(table(Clinical_data_ann[, IMS_group])), "BasalMyo_BasalHer2")

All_survival_analysis_data = data.frame(IMS = IMS_categories, p_value = 0, HR=0 , CI_lower=0, CI_upper = 0)
N.sets = length(IMS_categories)

dim(Clinical_data_ann)

if(TNBC == "TNBC"){
  Clinical_data_ann = Clinical_data_ann[which(Clinical_data_ann$bcr_patient_barcode %in% triple_negative$X),]
}

Survival_data_1 = Clinical_data_ann
Survival_data_1 = Survival_data_1[-which(is.na(Survival_data_1[, Ethnicity_variable])),]
if(Ethnicity_variable == "Ethnicity_reported_simple"){
  Survival_data_1 = Survival_data_1[-which(Survival_data_1[, Ethnicity_variable] %in% c("Asian", "Hispanic","Other")),]
}
if(Ethnicity_variable == "Assigned_Ethnicity_simplified"){
  Survival_data_1 = Survival_data_1[-which(Survival_data_1[, Ethnicity_variable] %in% c("Asian", "Unclear")),]
}


if(ICR_cluster %in% c("All", "ICR Medium-Low")){}else{
  Survival_data_1 = Survival_data_1[which(Survival_data_1$HML.ICR.Cluster == ICR_cluster),]
}

if(ICR_cluster == "ICR Medium-Low"){
  Survival_data_1$HML.ICR.Cluster[which(Survival_data_1$HML.ICR.Cluster %in% c("ICR Medium", "ICR Low"))] = "ICR Medium-Low"
  Survival_data_1 = Survival_data_1[which(Survival_data_1$HML.ICR.Cluster == "ICR Medium-Low"),]
}
dim(Survival_data_1)

if(Stages == "StageI_II"){
  Survival_data_1 = Survival_data_1[which(Survival_data_1$ajcc_pathologic_tumor_stage %in% c("Stage I", "Stage IA",
                                                                                             "Stage II", "Stage IIA", "Stage IIB")),]
}

if(Stages == "StageIII_IV"){
  Survival_data_1 = Survival_data_1[which(Survival_data_1$ajcc_pathologic_tumor_stage %in% c("Stage III", "Stage IIIA", "Stage IIIB",
                                                                                             "Stage IIIC", "Stage IV")),]
}

if(Stage_translation == "binary"){
  Survival_data_1$ajcc_pathologic_tumor_stage[which(Survival_data_1$ajcc_pathologic_tumor_stage %in% c("Stage I", "Stage IA",
                                                                                                       "Stage II", "Stage IIA", "Stage IIB"))] = "StageI_II"
  
  Survival_data_1$ajcc_pathologic_tumor_stage[which(Survival_data_1$ajcc_pathologic_tumor_stage %in% c("Stage III", "Stage IIIA", "Stage IIIB",
                                                                                                       "Stage IIIC", "Stage IV"))] = "StageIII_IV"
}

if(Stage_translation == "ordinal"){
  Survival_data_1$ajcc_pathologic_tumor_stage[which(Survival_data_1$ajcc_pathologic_tumor_stage %in% c("Stage I", "Stage IA"))] = 1
  Survival_data_1$ajcc_pathologic_tumor_stage[which(Survival_data_1$ajcc_pathologic_tumor_stage %in% c("Stage II", "Stage IIA", "Stage IIB"))] = 2
  Survival_data_1$ajcc_pathologic_tumor_stage[which(Survival_data_1$ajcc_pathologic_tumor_stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC"))] = 3
  Survival_data_1$ajcc_pathologic_tumor_stage[which(Survival_data_1$ajcc_pathologic_tumor_stage %in% c("Stage IV"))] = 4
}

#Survival_data_1 = Survival_data_1[which(Survival_data_1$bcr_patient_barcode %in% NeoAgs_Observed_Expected_patients),]
## Perform survival analysis

i = 2
for (i in 2:N.sets) {
  IMS = IMS_categories[i]
  colnames(Survival_data_1)[which(colnames(Survival_data_1) == "death_days_to")] = "days_to_death"
  colnames(Survival_data_1)[which(colnames(Survival_data_1) == "last_contact_days_to")] = "days_to_last_followup"
  
  if(IMS == "All"){
    Survival_data = Survival_data_1
  }else{
    if(IMS == "BasalMyo_BasalHer2"){
      Survival_data = Survival_data_1[which(Survival_data_1[,IMS_group] %in% c("BasalMyo", "BasalHer2")),]
    }else{
      Survival_data = Survival_data_1[which(Survival_data_1[, IMS_group] == IMS),]
    }
  }
  
  Y = Surv_cutoff_years * 365
  TS.EventFree = Survival_data[Survival_data[, Outcome] == "0", c(Outcome, paste0(Outcome, ".time"), Ethnicity_variable, 
                                                                  "ajcc_pathologic_tumor_stage", "age_at_initial_pathologic_diagnosis")]
  colnames(TS.EventFree) = c("Status","Time", "Group", "Stage", "Age")
  TS.EventFree$Time = as.numeric(as.character(TS.EventFree$Time))
  TS.EventFree$Time[TS.EventFree$Time > Y] = Y
  
  TS.EventOccured = Survival_data[Survival_data[, Outcome] == "1", c(Outcome, paste0(Outcome, ".time"), Ethnicity_variable, 
                                                                     "ajcc_pathologic_tumor_stage", "age_at_initial_pathologic_diagnosis")]
  colnames(TS.EventOccured) = c("Status","Time", "Group", "Stage", "Age")
  TS.EventOccured$Time = as.numeric(as.character(TS.EventOccured$Time))
  TS.EventOccured$Status[which(TS.EventOccured$Time> Y)] = "EventFree"
  TS.EventOccured$Time[TS.EventOccured$Time > Y] = Y
  
  TS.Surv = rbind (TS.EventOccured,TS.EventFree)
  TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
  TS.Surv$Status <- TS.Surv$Status == "1"
  TS.Surv = subset(TS.Surv,TS.Surv$Time > 1)                                                                                         # remove patients with less then 1 day follow up time
  
  # Number of evenets
  sum(TS.Surv$Status)
  
  # survival curve
  msurv = Surv(TS.Surv$Time/30.4, TS.Surv$Status)                                                                                    # calculate the number of months
  mfit = survfit(msurv~TS.Surv$Group,conf.type = "log-log")
  
  summary(mfit) 
  
  # Calculations
  mdiff = survdiff(eval(mfit$call$formula), data = eval(mfit$call$data))
  pval = pchisq(mdiff$chisq,length(mdiff$n) - 1,lower.tail = FALSE)
  pvaltxt = ifelse(pval < 0.0001,"p < 0.0001",paste("p =", signif(pval, 3)))
  
  #TS.Surv[,"Group"] = factor(TS.Surv[,"Group"], levels = c("ICR High", "ICR Medium", "ICR Low"))
  TS.Surv[,"Group"] = factor(TS.Surv[,"Group"], levels = c("White", "Black"))
  
  # Check this!!
  ##TS.Surv[,"Group"] = relevel(TS.Surv[,"Group"], "ICR High")
  mHR = coxph(formula = msurv ~ TS.Surv[,"Group"],data = TS.Surv, subset = TS.Surv$Group %in% c("White", "Black"))
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
  
  logrank_test(Surv(Time, Status) ~ Group, data = TS.Surv)
  summary(mHR)
  # plots
  #png(paste0("./Figures/Kaplan_meiers/C7_survival_by_IMS/Feb_2020_", TNBC, "_", Stages, "_", Ethnicity_variable,"_", ICR_cluster, "_", IMS_group, "_", IMS, "_", Outcome, "_Kaplan_Meijer_by_ethnicity.png"),
  #res=600, height=3,width=4,unit="in")                                                                                           # set filename
  # ggkm(mfit,
  #     timeby=12,
  #     ystratalabs = levels(TS.Surv[,"Group"]),
  #     ystrataname = NULL,
  #     main= paste0(IMS, "\n"),
  #     xlabs = "Time in months",
  #     cbPalette = cbPalette)
  #PLOT_HR = PLOT_HR,
  #PLOT_P = PLOT_P,
  #PLOT_CI1 = PLOT_CI1,
  #PLOT_CI2 = PLOT_CI2)
  #dev.off()
  plot = ggsurvplot(mfit,
                    data = TS.Surv,
                    censor = TRUE,
                    risk.table = TRUE,
                    tables.y.text.col = FALSE,
                    tables.y.text = FALSE,
                    tables.height = 0.3,
                    tables.theme = theme_cleantable(),
                    #tables.col = "strata",
                    risk.table.pos = "out",
                    legend = "none",
                    ylab = "",
                    xlab = "Time in months",
                    fontsize = 4.5,
                    font.x = 18,
                    font.tickslab = 18,
                    #censor.shape = 3,
                    #censor.size = 1,
                    #pval = TRUE
                    palette = cbPalette
  )
  png(paste0("./Figures/Kaplan_meiers/C7_survival_by_IMS/Aug_2020_", TNBC, "_", Stages, "_", Ethnicity_variable,"_", ICR_cluster, "_", IMS_group, "_", IMS, "_", Outcome, "_Kaplan_Meijer_by_ethnicity.png"),
      res=600, height=3.8,width=4.2,unit="in")                                                                                           # set filename
  print(plot)
  dev.off()     
  #PLOT_HR = PLOT_HR,
  #PLOT_P = PLOT_P,
  #PLOT_CI1 = PLOT_CI1,
  #PLOT_CI2 = PLOT_CI2)
  #####
}

#uni_variate = coxph(formula = msurv ~ TS.Surv[,"Group"],data = TS.Surv, subset = TS.Surv$Group %in% c("White", "Black"))

TS.Surv$Stage[which(TS.Surv$Stage %in% c("Stage X", "[Not Available]", "[Discrepancy]"))] = NA
TS.Surv$Stage = as.numeric(TS.Surv$Stage)

TS.Surv$Age = as.numeric(TS.Surv$Age)

TS.Surv$Stage = as.character(TS.Surv$Stage)
TS.Surv$Stage[which(TS.Surv$Stage %in% c("1", "2"))] = "StageI&II"
TS.Surv$Stage[which(TS.Surv$Stage %in% c("3", "4"))] = "StageIII&IV"

TS.Surv = TS.Surv[-which(is.na(TS.Surv$Stage)),]
multivariate_rev = coxph(formula = Surv(Time, Status) ~ Group + Stage + Age, data = TS.Surv)
cox.zph(multivariate_rev)
summary(multivariate_rev)

univariate = coxph(formula = Surv(Time, Status) ~ Group, data = TS.Surv)
summary(univariate)

univariate = coxph(formula = Surv(Time, Status) ~ Stage, data = TS.Surv)
summary(univariate)

univariate = coxph(formula = Surv(Time, Status) ~ Age, data = TS.Surv)
summary(univariate)
