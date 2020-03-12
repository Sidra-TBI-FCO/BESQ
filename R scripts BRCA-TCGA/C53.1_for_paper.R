
## Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA")                                                                   # Setwd to location were output files have to be saved.

source("../tools/ipak.function.R")
source("../tools/ggkm_C12.R")

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
Ethnicity = "White" # "Black"
only_ICR_Medium_and_Low = "only_ICR_Medium_and_Low"

# Load data
#load(paste0(code_path, "Datalists/ICR_genes.RData")) 
load("./Analysis/Sample_annotations.Rdata")
max_splits = read.csv(file = paste0("./Analysis/Hossam_prediction_model/Final_SHAP_p_value_and_cutoffs/", tolower(Ethnicity) ,"_bindea.csv"), stringsAsFactors = FALSE)
load("./Analysis/ssGSEA_scores/Bindea_ORIG_ES_scores.Rdata")
load("./Analysis/ICR data/TCGA_BRCA_table_cluster_assignment.RData")

# Create folders
dir.create("./Analysis/",showWarnings = FALSE)                                                                        # Create folder to save processed data (by Assembler module B)
dir.create("./Figures/Kaplan_meiers", showWarnings = FALSE)
dir.create("./Figures/Kaplan_meiers/C53.1_survival_by_pathway_category", showWarnings = FALSE)

max_splits$Feature = gsub("_", " ", max_splits$Feature)

rownames(ES) = gsub("\\/", " ", rownames(ES))
# Subset Survival data to only get BasalMyo and the above specified ethnicity
dim(Clinical_data_ann)
if(only_ICR_Medium_and_Low == "only_ICR_Medium_and_Low"){
  Clinical_data_ann = Clinical_data_ann[which(Clinical_data_ann$HML.ICR.Cluster %in% c("ICR Medium", "ICR Low")),]
}
Survival_data_1 = Clinical_data_ann
if(Ethnicity_variable == "Assigned_Ethnicity_simplified"){
  Survival_data_1 = Survival_data_1[-which(Survival_data_1[, Ethnicity_variable] %in% c("Asian", "Unclear")),]
}
dim(Survival_data_1)
colnames(Survival_data_1)[which(colnames(Survival_data_1) == "death_days_to")] = "days_to_death"
colnames(Survival_data_1)[which(colnames(Survival_data_1) == "last_contact_days_to")] = "days_to_last_followup"
if(IMS == "BasalMyo"){
  Survival_data = Survival_data_1[which(Survival_data_1[,IMS_group] == "BasalMyo"),]
}
if(Ethnicity == "All"){}else{
  Survival_data = Survival_data[which(Survival_data$Assigned_Ethnicity_simplified == Ethnicity),]
}

## Start loop (in loop the clusters are separated based on cutoff and survival plot is created accordingly)
## Survival data results are stored in dataframe "results"

pathways = max_splits$Feature
pathways = c("Th2 cells", "TReg", "DC", "B cells")
N.pathways = length(pathways)

results = data.frame(Pathway = pathways, HR = NA, P = NA, CI_Lower = NA, CI_Higher =NA)

i=3
for (i in 1:N.pathways){
  Pathway = pathways[i]
  plot_df = data.frame(Patient = Survival_data$bcr_patient_barcode,
                       OS = Survival_data$OS,
                       OS.Time = Survival_data$OS.time,
                       Pathway_score = NA,
                       Pathway_category = NA)
  if(Pathway == "ICRscore"){
    plot_df$Pathway_score = clustering$ICRscore[match(plot_df$Patient, substring(rownames(clustering), 1, 12))]
  }else{
    plot_df$Pathway_score = ES[Pathway,][match(plot_df$Patient, substring(colnames(ES), 1, 12))]
  }
  #cutoff = max_splits$Bound[which(max_splits$Feature == Pathway)]
  cutoff = median(plot_df$Pathway_score)
  plot_df$Pathway_category[which(plot_df$Pathway_score < cutoff)] = "Pathway Low"
  plot_df$Pathway_category[which(plot_df$Pathway_score >= cutoff)] = "Pathway High"
  plot_df$OS = as.character(plot_df$OS)
  
  Y = Surv_cutoff_years * 365
  TS.EventFree = plot_df[which(plot_df$OS == "0"), c("OS", "OS.Time", "Pathway_category")]
  colnames(TS.EventFree) = c("Status","Time", "Group")
  TS.EventFree$Time = as.numeric(as.character(TS.EventFree$Time))
  TS.EventFree$Time[TS.EventFree$Time > Y] = Y
  
  TS.EventOccured = plot_df[which(plot_df$OS == "1"), c("OS", "OS.Time", "Pathway_category")]
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
  if(classification == "binary"){
    TS.Surv[,"Group"] = factor(TS.Surv[,"Group"], levels = c("Pathway High", "Pathway Low"))
  }
  
  # Check this!!
  ##TS.Surv[,"Group"] = relevel(TS.Surv[,"Group"], "ICR Pathway High")
  mHR = coxph(formula = msurv ~ TS.Surv[,"Group"],data = TS.Surv, subset = TS.Surv$Group %in% c("Pathway High", "Pathway Low"))
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
  if(classification == "binary"){
    PLOT_P = signif(p[1], digits = 3)
    PLOT_HR = round(signif(exp(mHR.extract@coef),3)[1], 3)
    PLOT_CI1 = CI[1,1]
    PLOT_CI2 = CI[1,2]
  }
  
  results[which(results$Pathway == Pathway),2:5] = c(PLOT_HR, PLOT_P, PLOT_CI1, PLOT_CI2)
  
  # plots
  png(paste0("./Figures/Kaplan_meiers/C53.1_survival_by_pathway_category/", only_ICR_Medium_and_Low, "_Median_Oct_2019_", Ethnicity_variable,"_", Ethnicity, "_", Pathway,"_", classification, "_", Outcome, "_Kaplan_Meijer_by_cluster.png"),
      res=600, height=3,width=6,unit="in")                                                                                           # set filename
  ggkm(mfit,
       timeby=12,
       ystratalabs = levels(TS.Surv[,"Group"]),
       ystrataname = NULL,
       main= paste0(IMS, " ", Pathway, " with cutoff ", round(cutoff, 5) ," in ", Ethnicity),
       xlabs = "Time in months",
       Title_size = 18,
       palette = c("#FF00FF", "#40E0D0"),
       legend = "none")
       #PLOT_HR = PLOT_HR,
       #PLOT_P = PLOT_P,
       #PLOT_CI1 = PLOT_CI1,
       #PLOT_CI2 = PLOT_CI2)
  dev.off()
  
  
}
#dir.create("Analysis/C51_results_ggkm_optimal_cutoffs", showWarnings = FALSE)
#save(results, file = paste0("./Analysis/C51_results_ggkm_optimal_cutoffs/C51_", Ethnicity, "_", IMS, "_results_ggkm.Rdata"))
#write.csv(results, file = paste0("./Analysis/C51_results_ggkm_optimal_cutoffs/C51_", Ethnicity, "_", IMS, "_results_ggkm.csv"))






