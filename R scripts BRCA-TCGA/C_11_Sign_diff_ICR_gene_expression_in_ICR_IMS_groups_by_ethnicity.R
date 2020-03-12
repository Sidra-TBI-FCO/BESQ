

# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA")                                                                   # Setwd to location were output files have to be saved.

source("../tools/ipak.function.R")

ipak(c("ggplot2", "ggpubr", "dplyr"))

# Set parameters
Ethnicity_variable = "Assigned_Ethnicity_simplified" #"Assigned_Ethnicity_simplified" #"Ethnicity_reported_simple"

# Load data
load("./Analysis/Sample_annotations.Rdata")
load("./Data/expression matrix/BRCA_gene_RNAseq_normalized_Tissue_Type_and_whitelists_filtered.Rdata")
load("../tools/ICR_genes.RData")

dir.create("./Figures/C11_Boxplot_ICR_gene_across_ethnicities", showWarnings = FALSE)

# Analysis
RNAseq_log2 = log(filtered.norm.RNAseqData +1, 2)

# Delete NA
Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann[, Ethnicity_variable])),]
if(Ethnicity_variable == "Ethnicity_reported_simple"){
  Clinical_data_ann = Clinical_data_ann[-which(Clinical_data_ann[, Ethnicity_variable] == "Other"),]
}
if(Ethnicity_variable == "Assigned_Ethnicity_simplified"){
  Clinical_data_ann = Clinical_data_ann[-which(Clinical_data_ann[, Ethnicity_variable] == "Unclear"),]
}

Clinical_data_ann = Clinical_data_ann[-which(is.na(Clinical_data_ann$IMS)),]

ICR_clusters = names(table(Clinical_data_ann$HML.ICR.Cluster))
IMS_categories = names(table(Clinical_data_ann$IMS_Mathews))
Ethicities = names(table(Clinical_data_ann[, Ethnicity_variable]))

p_value_df = Clinical_data_ann %>%
  group_by(HML.ICR.Cluster, IMS_Mathews) %>%
  summarise(count=n())

p_value_df$count = NULL
p_value_df$gene = NA

p_value_df_final = p_value_df
p_value_df_final$gene = ICR_genes[1]

for (i in 2:length(ICR_genes)){
  gene = ICR_genes[i]
  p_value_df$gene = gene
  p_value_df_final = rbind(p_value_df_final, p_value_df)
}

p_value_df_final[,c("Mean White", "Mean Black", "Mean Asian",
                    "p_value_W_B", "p_value_W_A", "p_value_B_A")] = NA

i = 1
j = 2
k = 1
for (i in 1:length(ICR_clusters)){
  ICR_cluster = ICR_clusters[i]
  for (j in 1:length(IMS_categories)){
    IMS = IMS_categories[j]
    for (k in 1:length(ICR_genes)){
      gene = ICR_genes[k]
      
      plot_df = Clinical_data_ann[, c("bcr_patient_barcode","HML.ICR.Cluster", "IMS_Mathews")]
      plot_df$Ethnicity = Clinical_data_ann[, Ethnicity_variable]
      plot_df$gene = RNAseq_log2[gene,][match(plot_df$bcr_patient_barcode,
                                              substring(colnames(RNAseq_log2), 1, 12))]
      plot_df = plot_df[which(plot_df$HML.ICR.Cluster == ICR_cluster & plot_df$IMS_Mathews == IMS),]
      
      Ethnicities_with_small_groups = names(which(table(plot_df$Ethnicity) < 3))
      
      agg = aggregate(plot_df$gene,
                      by = list(plot_df$Ethnicity),
                      FUN = mean)
      means_df = data.frame(Ethnicity = c("White", "Black", "Asian"), Mean = NA)
      means_df$Mean = agg$x[match(means_df$Ethnicity, agg$Group.1)]
      
      if(sum(c("White", "Black") %in% Ethnicities_with_small_groups) > 0){
        p1 = NA
      }else{ p1 = round(t.test(plot_df[which(plot_df$Ethnicity %in% c("White", "Black")),]$gene ~ plot_df[which(plot_df$Ethnicity %in% c("White", "Black")),]$Ethnicity,
                               data = plot_df[which(plot_df$Ethnicity %in% c("White", "Black")),], paired = FALSE, var.equal = FALSE)$p.value, 3)
      }
      if(sum(c("White", "Asian") %in% Ethnicities_with_small_groups) > 0){
        p2 = NA
      }else{p2 = round(t.test(plot_df[which(plot_df$Ethnicity %in% c("White", "Asian")),]$gene ~ plot_df[which(plot_df$Ethnicity %in% c("White", "Asian")),]$Ethnicity,
                              data = plot_df[which(plot_df$Ethnicity %in% c("White", "Asian")),], paired = FALSE, var.equal = FALSE)$p.value, 3)
      }
     # if(sum(c("White", "Hispanic") %in% Ethnicities_with_small_groups) > 0){
      #  p3 = NA
     # }else{p3 = round(t.test(plot_df[which(plot_df$Ethnicity %in% c("White", "Hispanic")),]$gene ~ plot_df[which(plot_df$Ethnicity %in% c("White", "Hispanic")),]$Ethnicity,
      #                        data = plot_df[which(plot_df$Ethnicity %in% c("White", "Hispanic")),], paired = FALSE, var.equal = FALSE)$p.value, 3)}
      if(sum(c("Black", "Asian") %in% Ethnicities_with_small_groups) > 0){
        p4 = NA
      }else{p4 = round(t.test(plot_df[which(plot_df$Ethnicity %in% c("Black", "Asian")),]$gene ~ plot_df[which(plot_df$Ethnicity %in% c("Black", "Asian")),]$Ethnicity,
                              data = plot_df[which(plot_df$Ethnicity %in% c("Black", "Asian")),], paired = FALSE, var.equal = FALSE)$p.value, 3)}
      #if(sum(c("Black", "Hispanic") %in% Ethnicities_with_small_groups) > 0){
      #  p5 = NA
      #}else{p5 = round(t.test(plot_df[which(plot_df$Ethnicity %in% c("Black", "Hispanic")),]$gene ~ plot_df[which(plot_df$Ethnicity %in% c("Black", "Hispanic")),]$Ethnicity,
      #                        data = plot_df[which(plot_df$Ethnicity %in% c("Black", "Hispanic")),], paired = FALSE, var.equal = FALSE)$p.value, 3)}
      #if(sum(c("Asian", "Hispanic") %in% Ethnicities_with_small_groups) > 0){
       # p6 = NA
      #}else{p6 = round(t.test(plot_df[which(plot_df$Ethnicity %in% c("Asian", "Hispanic")),]$gene ~ plot_df[which(plot_df$Ethnicity %in% c("Asian", "Hispanic")),]$Ethnicity,
     #                         data = plot_df[which(plot_df$Ethnicity %in% c("Asian", "Hispanic")),], paired = FALSE, var.equal = FALSE)$p.value, 3)}
      
      
      
      p_value_df_final[which(p_value_df_final$HML.ICR.Cluster == ICR_cluster & p_value_df_final$IMS_Mathews == IMS & p_value_df_final$gene == gene),
                       c("p_value_W_B", "p_value_W_A", "p_value_B_A")] = c(p1, p2, p4)
      p_value_df_final[which(p_value_df_final$HML.ICR.Cluster == ICR_cluster & p_value_df_final$IMS_Mathews == IMS & p_value_df_final$gene == gene),
                       c("Mean White", "Mean Black", "Mean Asian")] = means_df$Mean
      
      
      my_comparisons = list(c("White", "Black"), c("White", "Asian"),
                            c("Black", "Asian"))
      
      plot = ggplot(plot_df, aes(x = Ethnicity, y = gene)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.3, size = 0.5) +
        stat_compare_means(comparisons = my_comparisons, method = "t.test") +
        ylab(gene) +
        theme_bw() +
        ggtitle(paste0(gene, " expression across ethnicities\n",
                       "in ", IMS, " breast cancer in ", ICR_cluster, " cluster"))
      
      #png(filename = paste0("./Figures/C11_Boxplot_ICR_gene_across_ethnicities/Boxplot_", ICR_cluster,
                        #    "_", IMS, "_", gene, "_by_ethnicities.png"), width = 5, height = 5,
          #units = "in", res = 600)
     # plot(plot)
     # dev.off()
    }
  }
}

dir.create("./Analysis/t_test_between_ethnicities", showWarnings = FALSE)

# Analysis on p-values and directionality
p_value_df_final$Black_vs_White = p_value_df_final$`Mean White` > p_value_df_final$`Mean Black`
p_value_df_final$Black_vs_White[which(p_value_df_final$Black_vs_White == TRUE)] = "Down"
p_value_df_final$Black_vs_White[which(p_value_df_final$Black_vs_White == FALSE)] = "Up"

dir.create("./Analysis/t_test_between_ethnicities", showWarnings = FALSE)

signif = p_value_df_final[which(p_value_df_final$p_value_W_B < 0.05),]
vars_of_interest = table(signif$gene)[which(table(signif$gene) > 0)]
vars_of_interest = data.frame(vars_of_interest)
vars_of_interest$N_up_sig = NA
vars_of_interest$N_down_sig = NA
vars_of_interest$N_up = NA
vars_of_interest$N_down = NA


i=1
for (i in 1:nrow(vars_of_interest)){
  gene = vars_of_interest$Var1[i]
  subset = p_value_df_final[which(p_value_df_final$gene == gene),]
  sig = subset[which(subset$p_value_W_B < 0.05),]
  N_up_sig = sum(sig$Black_vs_White == "Up")
  N_down_sig = sum(sig$Black_vs_White == "Down")
  N_up = sum(subset$Black_vs_White == "Up", na.rm = TRUE)
  N_down = sum(subset$Black_vs_White == "Down", na.rm = TRUE)
  vars_of_interest$N_up_sig[which(vars_of_interest$Var1 == gene)] = N_up_sig
  vars_of_interest$N_down_sig[which(vars_of_interest$Var1 == gene)] = N_down_sig
  vars_of_interest$N_up[which(vars_of_interest$Var1 == gene)] = N_up
  vars_of_interest$N_down[which(vars_of_interest$Var1 == gene)] = N_down
}

vars_of_interest = vars_of_interest[order(vars_of_interest$Freq, decreasing = TRUE),]
vars_of_interest$ratio = vars_of_interest$N_up / vars_of_interest$N_down
vars_of_interest = vars_of_interest[order(vars_of_interest$ratio, decreasing = TRUE),]

write.csv(vars_of_interest, file = paste0("./Analysis/t_test_between_ethnicities/C11_", Ethnicity_variable, "_Diff_", "ICR_genes","_between_ICR_IMS_Ethnicity.csv"))
save(p_value_df_final, vars_of_interest, file = paste0("./Analysis/t_test_between_ethnicities/C11_", Ethnicity_variable, "_Diff_", "ICR_genes", "_between_ICR_IMS_Ethnicity.Rdata"))
