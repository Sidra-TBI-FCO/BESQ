
### SCRIPT NOT SUITABLE FOR RA-QA Cohort


# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-RA-QA")                                                                   # Setwd to location were output files have to be saved.

source("../tools/ipak.function.R")

ipak(c("ggplot2", "ggpubr", "dplyr"))

# Set parameters
Gene_set = "Bindea_ORIG"

# Load data
load("./Analysis/Sample_annotations.Rdata")
load(paste0("./Analysis/ssGSEA_scores/", Gene_set,"_ES_scores.Rdata"))
#load("./Analysis/ssGSEA_scores/Selected.pathways_ES_scores.Rdata")

# Analysis
cell_types = rownames(ES)

# Delete NA
Clinical_data_ann = Clinical_data_ann[-which(Clinical_data_ann$Ethnicity == "Other"),]
Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann$IMS)),]

ICR_clusters = names(table(Clinical_data_ann$HL.ICR.Cluster))
IMS_categories = names(table(Clinical_data_ann$IMS))
Ethicities = names(table(Clinical_data_ann$Ethnicity))

p_value_df = Clinical_data_ann %>%
  group_by(HL.ICR.Cluster, IMS) %>%
  summarise(count=n())

p_value_df$count = NULL
p_value_df$cell_type = NA

p_value_df_final = p_value_df
p_value_df_final$cell_type = cell_types[1]

for (i in 2:length(cell_types)){
  cell_type = cell_types[i]
  p_value_df$cell_type = cell_type
  p_value_df_final = rbind(p_value_df_final, p_value_df)
}

p_value_df_final[,c("Mean White", "Mean Black", "Mean Asian", "Mean Hispanic", 
                    "p_value_W_B", "p_value_W_A", "p_value_W_H", "p_value_B_A", "p_value_B_H", "p_value_A_H")] = NA

i = 1
j = 1
k = 1
for (i in 1:length(ICR_clusters)){
  ICR_cluster = ICR_clusters[i]
  for (j in 1:length(IMS_categories)){
    IMS = IMS_categories[j]
    for (k in 1:length(cell_types)){
      cell_type = cell_types[k]
      
      plot_df = Clinical_data_ann[, c("Sample_ID","HL.ICR.Cluster", "IMS", "Ethnicity")]
      plot_df$cell_type = ES[cell_type,][match(plot_df$Sample_ID,
                                               colnames(ES))]
      plot_df = plot_df[which(plot_df$HL.ICR.Cluster == ICR_cluster & plot_df$IMS == IMS),]
      
      Ethnicities_with_small_groups = names(which(table(plot_df$Ethnicity) < 3))
      
      agg = aggregate(plot_df$cell_type,
                      by = list(plot_df$Ethnicity),
                      FUN = mean)
      means_df = data.frame(Ethnicity = c("White", "Black", "Asian", "Hispanic"), Mean = NA)
      means_df$Mean = agg$x[match(means_df$Ethnicity, agg$Group.1)]
      
      if(sum(c("White", "Black") %in% Ethnicities_with_small_groups) > 0){
        p1 = NA
      }else{ p1 = round(t.test(plot_df[which(plot_df$Ethnicity %in% c("White", "Black")),]$cell_type ~ plot_df[which(plot_df$Ethnicity %in% c("White", "Black")),]$Ethnicity,
                               data = plot_df[which(plot_df$Ethnicity %in% c("White", "Black")),], paired = FALSE, var.equal = FALSE)$p.value, 3)
      }
      if(sum(c("White", "Asian") %in% Ethnicities_with_small_groups) > 0){
        p2 = NA
      }else{p2 = round(t.test(plot_df[which(plot_df$Ethnicity %in% c("White", "Asian")),]$cell_type ~ plot_df[which(plot_df$Ethnicity %in% c("White", "Asian")),]$Ethnicity,
                              data = plot_df[which(plot_df$Ethnicity %in% c("White", "Asian")),], paired = FALSE, var.equal = FALSE)$p.value, 3)
      }
      if(sum(c("White", "Hispanic") %in% Ethnicities_with_small_groups) > 0){
        p3 = NA
      }else{p3 = round(t.test(plot_df[which(plot_df$Ethnicity %in% c("White", "Hispanic")),]$cell_type ~ plot_df[which(plot_df$Ethnicity %in% c("White", "Hispanic")),]$Ethnicity,
                              data = plot_df[which(plot_df$Ethnicity %in% c("White", "Hispanic")),], paired = FALSE, var.equal = FALSE)$p.value, 3)}
      if(sum(c("Black", "Asian") %in% Ethnicities_with_small_groups) > 0){
        p4 = NA
      }else{p4 = round(t.test(plot_df[which(plot_df$Ethnicity %in% c("Black", "Asian")),]$cell_type ~ plot_df[which(plot_df$Ethnicity %in% c("Black", "Asian")),]$Ethnicity,
                              data = plot_df[which(plot_df$Ethnicity %in% c("Black", "Asian")),], paired = FALSE, var.equal = FALSE)$p.value, 3)}
      if(sum(c("Black", "Hispanic") %in% Ethnicities_with_small_groups) > 0){
        p5 = NA
      }else{p5 = round(t.test(plot_df[which(plot_df$Ethnicity %in% c("Black", "Hispanic")),]$cell_type ~ plot_df[which(plot_df$Ethnicity %in% c("Black", "Hispanic")),]$Ethnicity,
                              data = plot_df[which(plot_df$Ethnicity %in% c("Black", "Hispanic")),], paired = FALSE, var.equal = FALSE)$p.value, 3)}
      if(sum(c("Asian", "Hispanic") %in% Ethnicities_with_small_groups) > 0){
        p6 = NA
      }else{p6 = round(t.test(plot_df[which(plot_df$Ethnicity %in% c("Asian", "Hispanic")),]$cell_type ~ plot_df[which(plot_df$Ethnicity %in% c("Asian", "Hispanic")),]$Ethnicity,
                              data = plot_df[which(plot_df$Ethnicity %in% c("Asian", "Hispanic")),], paired = FALSE, var.equal = FALSE)$p.value, 3)}
      
      
      
      p_value_df_final[which(p_value_df_final$HML.ICR.Cluster == ICR_cluster & p_value_df_final$IMS == IMS & p_value_df_final$cell_type == cell_type),
                       c("p_value_W_B", "p_value_W_A", "p_value_W_H", "p_value_B_A", "p_value_B_H", "p_value_A_H")] = c(p1, p2, p3, p4, p5, p6)
      p_value_df_final[which(p_value_df_final$HML.ICR.Cluster == ICR_cluster & p_value_df_final$IMS == IMS & p_value_df_final$cell_type == cell_type),
                       c("Mean White", "Mean Black", "Mean Asian", "Mean Hispanic")] = means_df$Mean
      
    }
  }
}


# Analysis on p-values and directionality
p_value_df_final$Black_vs_White = p_value_df_final$`Mean White` > p_value_df_final$`Mean Black`
p_value_df_final$Black_vs_White[which(p_value_df_final$Black_vs_White == TRUE)] = "Down"
p_value_df_final$Black_vs_White[which(p_value_df_final$Black_vs_White == FALSE)] = "Up"

dir.create("./Analysis/t_test_between_ethnicities", showWarnings = FALSE)

signif = p_value_df_final[which(p_value_df_final$p_value_W_B < 0.05),]
vars_of_interest = table(signif$cell_type)[which(table(signif$cell_type) > 0)]
vars_of_interest = data.frame(vars_of_interest)
vars_of_interest$N_up_sig = NA
vars_of_interest$N_down_sig = NA
vars_of_interest$N_up = NA
vars_of_interest$N_down = NA


i=1
for (i in 1:nrow(vars_of_interest)){
  var = vars_of_interest$Var1[i]
  subset = p_value_df_final[which(p_value_df_final$cell_type == var),]
  sig = subset[which(subset$p_value_W_B < 0.05),]
  N_up_sig = sum(sig$Black_vs_White == "Up")
  N_down_sig = sum(sig$Black_vs_White == "Down")
  N_up = sum(subset$Black_vs_White == "Up")
  N_down = sum(subset$Black_vs_White == "Down")
  vars_of_interest$N_up_sig[which(vars_of_interest$Var1 == var)] = N_up_sig
  vars_of_interest$N_down_sig[which(vars_of_interest$Var1 == var)] = N_down_sig
  vars_of_interest$N_up[which(vars_of_interest$Var1 == var)] = N_up
  vars_of_interest$N_down[which(vars_of_interest$Var1 == var)] = N_down
}

vars_of_interest = vars_of_interest[order(vars_of_interest$Freq, decreasing = TRUE),]
write.csv(vars_of_interest, file = paste0("./Analysis/t_test_between_ethnicities/C14_Diff_", Gene_set,"_between_ICR_IMS_Ethnicity.csv"))
save(p_value_df_final, vars_of_interest, file = paste0("./Analysis/t_test_between_ethnicities/C14_Diff_", Gene_set, "_between_ICR_IMS_Ethnicity.Rdata"))
