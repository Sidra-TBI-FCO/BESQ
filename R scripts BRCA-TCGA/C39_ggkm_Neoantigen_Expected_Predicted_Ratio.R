
## Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA")                                                                   # Setwd to location were output files have to be saved.

source("../tools/ipak.function.R")
source("../tools/ggkm_July_2019.R")

required.packages = c("survival","reshape","ggplot2","plyr","Rcpp","colorspace","texreg", "Hmisc")
required.bioconductor.packages = "survival"
ipak(required.packages)
ibiopak(required.bioconductor.packages)

# Set Parameters                                                                                            # "HML_classification" or "k3" or "k4" or "k5"
Surv_cutoff_years = 10
Outcome = "OS" # "DSS" "OS"
Ethnicity_variable = "Assigned_Ethnicity_simplified" #"Ethnicity_reported_simple"
cbPalette = c("lightblue", "orange")
IMS_group = "IMS_Mathews" # "IMS", "IMS_Mathews"
IMS = "BasalMyo"
classification = "binary" #"tertiles" # "binary"
Ethnicity = "All" # "Black"
#ICR_groups = c("ICR High")
cutoff = 1

# Load data
#load(paste0(code_path, "Datalists/ICR_genes.RData")) 
load("./Analysis/Sample_annotations.Rdata")
Rooney = read.csv("./Data/Rooney_expected_observed_neontigen_ratio-mmc4.csv",
                  stringsAsFactors = FALSE)

# Create folders
dir.create("./Analysis/",showWarnings = FALSE)                                                                        # Create folder to save processed data (by Assembler module B)
dir.create("./Figures/Kaplan_meiers", showWarnings = FALSE)
dir.create("./Figures/Kaplan_meiers/C39_survival_by_neoantigen_mutational_load_ratio", showWarnings = FALSE)

# Subset Survival data to only get BasalMyo and the above specified ethnicity
dim(Clinical_data_ann)
Survival_data_1 = Clinical_data_ann
if(Ethnicity_variable == "Assigned_Ethnicity_simplified"){
  Survival_data_1 = Survival_data_1[-which(Survival_data_1[, Ethnicity_variable] %in% c("Asian", "Unclear")),]
}
dim(Survival_data_1)
colnames(Survival_data_1)[which(colnames(Survival_data_1) == "death_days_to")] = "days_to_death"
colnames(Survival_data_1)[which(colnames(Survival_data_1) == "last_contact_days_to")] = "days_to_last_followup"
if(IMS == "BasalMyo"){
  Survival_data = Survival_data_1[which(Survival_data_1[,IMS_group] == "BasalMyo"),]
}else{Survival_data = Survival_data_1}
if(Ethnicity == "All"){}else{
  Survival_data = Survival_data[which(Survival_data$Assigned_Ethnicity_simplified == Ethnicity),]
}
Survival_data = Survival_data[which(Survival_data$HML.ICR.Cluster %in% ICR_groups),]

plot_df = data.frame(Patient = Survival_data$bcr_patient_barcode,
                     OS = Survival_data$OS,
                     OS.Time = Survival_data$OS.time,
                     NeoAgs_Observed_Expected = NA,
                     NeoAgs_Observed_Expected_category = NA)

plot_df$NeoAgs_Observed_Expected = Rooney$NeoAgs_Observed.Expected[match(plot_df$Patient,
                                                                         Rooney$PatientID)]
plot_df = plot_df[-which(is.na(plot_df$NeoAgs_Observed_Expected)),]

plot_df$NeoAgs_Observed_Expected_category[which(plot_df$NeoAgs_Observed_Expected < cutoff)] = "NeoAgs Observed/Expected Low"
plot_df$NeoAgs_Observed_Expected_category[which(plot_df$NeoAgs_Observed_Expected >= cutoff)] = "NeoAgs Observed/Expected High"
plot_df$OS = as.character(plot_df$OS)

NeoAgs_Observed_Expected_patients = as.character(plot_df$Patient)
save(NeoAgs_Observed_Expected_patients, file = "./Figures/Kaplan_meiers/C39_survival_by_neoantigen_mutational_load_ratio/NeoAgs_Observed_Expected_patients.Rdata")

Y = Surv_cutoff_years * 365
TS.EventFree = plot_df[which(plot_df$OS == "0"), c("OS", "OS.Time", "NeoAgs_Observed_Expected_category")]
colnames(TS.EventFree) = c("Status","Time", "Group")
TS.EventFree$Time = as.numeric(as.character(TS.EventFree$Time))
TS.EventFree$Time[TS.EventFree$Time > Y] = Y

TS.EventOccured = plot_df[which(plot_df$OS == "1"), c("OS", "OS.Time", "NeoAgs_Observed_Expected_category")]
colnames(TS.EventOccured) = c("Status","Time", "Group")
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

#TS.Surv[,"Group"] = factor(TS.Surv[,"Group"], levels = c("ICR Pathway High", "ICR Medium", "ICR Low"))

TS.Surv[,"Group"] = factor(TS.Surv[,"Group"], levels = c("NeoAgs Observed/Expected High", "NeoAgs Observed/Expected Low"))

# Check this!!
##TS.Surv[,"Group"] = relevel(TS.Surv[,"Group"], "ICR Pathway High")
mHR = coxph(formula = msurv ~ TS.Surv[,"Group"],data = TS.Surv, subset = TS.Surv$Group %in% c("NeoAgs Observed/Expected High", "NeoAgs Observed/Expected Low"))
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

## extract only ICR High vs. ICR Pathway Low

PLOT_P = signif(p[1], digits = 3)
PLOT_HR = round(signif(exp(mHR.extract@coef),3)[1], 3)
PLOT_CI1 = CI[1,1]
PLOT_CI2 = CI[1,2]  

# plots
png(paste0("./Figures/Kaplan_meiers/C39_survival_by_neoantigen_mutational_load_ratio/Kaplan_Meijer_by_Immunoediting_score_in_", IMS ,
           "_cutoff_", cutoff, "_cluster.png"),
    res=600, height=6,width=8,unit="in")                                                                                           # set filename
ggkm(mfit,
     timeby=12,
     ystratalabs = levels(TS.Surv[,"Group"]),
     ystrataname = NULL,
     main= paste0(IMS),
     xlabs = "Time in months",
     Title_size = 18,
     palette = c("purple", "orange"),
     legend = "none",
     PLOT_HR = PLOT_HR,
     PLOT_P = PLOT_P,
     PLOT_CI1 = PLOT_CI1,
     PLOT_CI2 = PLOT_CI2)

dev.off()




