
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
IMS = "All"
#classification = "binary" #"tertiles" # "binary"
Ethnicity = "All" # "Black"
ICR_groups = c("ICR High")
#cutoff = 1

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

Survival_data_1 = Survival_data_1[which(Survival_data_1$HML.ICR.Cluster %in% ICR_groups),]

dim(Survival_data_1)
colnames(Survival_data_1)[which(colnames(Survival_data_1) == "death_days_to")] = "days_to_death"
colnames(Survival_data_1)[which(colnames(Survival_data_1) == "last_contact_days_to")] = "days_to_last_followup"

subtypes = c("All", unique(as.character(Survival_data_1$IMS_Mathews)))
#subtypes = subtypes[-which(is.na(subtypes))]

results = data.frame(IMS = subtypes, HR = NA, p_value = NA, 
                     CI_lower = NA, CI_upper = NA,
                     N = NA)

i=1  
for (i in 1:length(subtypes)){
  IMS = subtypes[i]
  if(IMS == "All"){Survival_data = Survival_data_1}else{
    Survival_data = Survival_data_1[which(Survival_data_1[,IMS_group] == IMS),]
  }
  if(Ethnicity == "All"){}else{
    Survival_data = Survival_data[which(Survival_data$Assigned_Ethnicity_simplified == Ethnicity),]
  }
  
  plot_df = data.frame(Patient = Survival_data$bcr_patient_barcode,
                       OS = Survival_data$OS,
                       OS.Time = Survival_data$OS.time,
                       NeoAgs_Observed_Expected = NA)
  
  plot_df$NeoAgs_Observed_Expected = Rooney$NeoAgs_Observed.Expected[match(plot_df$Patient,
                                                                           Rooney$PatientID)]
  plot_df = plot_df[-which(is.na(plot_df$NeoAgs_Observed_Expected)),]
  
  plot_df$OS = as.character(plot_df$OS)
  
  Y = Surv_cutoff_years * 365
  TS.EventFree = plot_df[which(plot_df$OS == "0"), c("OS", "OS.Time", "NeoAgs_Observed_Expected")]
  colnames(TS.EventFree) = c("Status","Time", "Group")
  TS.EventFree$Time = as.numeric(as.character(TS.EventFree$Time))
  TS.EventFree$Time[TS.EventFree$Time > Y] = Y
  
  TS.EventOccured = plot_df[which(plot_df$OS == "1"), c("OS", "OS.Time", "NeoAgs_Observed_Expected")]
  colnames(TS.EventOccured) = c("Status","Time", "Group")
  TS.EventOccured$Time = as.numeric(as.character(TS.EventOccured$Time))
  TS.EventOccured$Status[which(TS.EventOccured$Time> Y)] = "0"
  TS.EventOccured$Time[TS.EventOccured$Time > Y] = Y
  
  TS.Surv = rbind (TS.EventOccured,TS.EventFree)
  TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
  TS.Surv$Status <- TS.Surv$Status == "1"
  TS.Surv = subset(TS.Surv,TS.Surv$Time > 1)  
  
  if(nrow(TS.EventFree) == 0 | nrow(TS.EventOccured) == 0){next}
  
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
  
  results$p_value[which(results$IMS == IMS)] = p
  results$HR[which(results$IMS == IMS)] = summary$coefficients[2]
  results$CI_lower[which(results$IMS == IMS)] = summary$conf.int[3]
  results$CI_upper[which(results$IMS == IMS)] = summary$conf.int[4]
  results$N[which(results$IMS == IMS)] = nrow(TS.Surv)
}

save(results, file = paste0("./Analysis/C39_coxph_Neoantigen_Expected_Observed_Ratio_", 
                            Ethnicity, "_", ICR_groups,".Rdata"))
