
# Setup environment
rm(list = ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA")  

source("../tools/ipak.function.R")
required.packages <- c("RColorBrewer", "forestplot")
ipak(required.packages)

# Set parameters
Outcome = "OS"
Ethnicity = "All"  # "All", "White", "Black"
IMS = "BasalMyo" # "BasalMyo"

# Load
HR.table = data.frame(Cell_types = c("All", "Stage I and II", "Stage III and IV"),
                      p_value = c(0.020,0.0905, 0.373),
                      HR = c(2.39,2.36, 1.66),
                      CI_lower = c(1.15,0.85,0.54),
                      CI_upper = c(4.98,6.55,5.09))

n.cells = nrow(HR.table)
x = n.cells + 2

HR.matrix = as.matrix(HR.table)
rownames(HR.matrix) = HR.matrix[,1]
HR.matrix = HR.matrix[,-c(1)]
mode(HR.matrix) = "numeric"

HR.table$p_value = signif(HR.table$p_value, 2)
HR.table$HR = signif(HR.table$HR, 3)
tabletext<-cbind(
  c("Cell", as.character(HR.table$Cell_types)[c(1:n.cells)]),
  c("p-value", HR.table$p_value[c(1:n.cells)]),
  c("HR",      HR.table$HR[c(1:n.cells)]))


dir.create("./Figures/C41_Forest_plot_by_stage", showWarnings = FALSE)
pdf(file = paste0("./Figures/C41_Forest_plot_by_stage/C41_stratified_by_stage",
                  IMS, "_", Ethnicity,".pdf"),
    height = 2, width = 3)

forestplot(mean = HR.matrix[,"HR"],
           lower = HR.matrix[,"CI_lower"],
           upper = HR.matrix[,"CI_upper"],
           labeltext = tabletext[-1,],
           new_page = FALSE,
           zero = 1,
           #is.summary=c(TRUE,rep(FALSE,n.cells),TRUE,rep(FALSE,n.cells),TRUE,FALSE),
           clip=c(0.001,55),
           xlog=TRUE,
           #xlim = c(0, 4),
           boxsize = .1,
           vertices = FALSE,
           col=fpColors(box="darkblue",line="darkgrey"),
           txt_gp = fpTxtGp(label = gpar(fontsize = 8.5), xlab = gpar(fontsize = 8.5),
                            ticks = gpar(fontsize = 10))
)
dev.off()


