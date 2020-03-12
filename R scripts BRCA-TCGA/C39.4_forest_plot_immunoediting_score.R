
# Setup environment
rm(list = ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA")  

source("../tools/ipak.function.R")
required.packages <- c("RColorBrewer", "forestplot")
ipak(required.packages)

# Set parameters
Outcome = "OS"
Ethnicity = "Black"  # "All", "White", "Black"

# Load
load(paste0("./Analysis/C39_coxph_Neoantigen_Expected_Observed_Ratio_", Ethnicity, ".Rdata"))

HR.table = results
#HR.table = HR.table[-which(HR.table$CI_lower == 0),]

n.cells = nrow(HR.table)
x = n.cells + 2


HR.matrix = as.matrix(HR.table)
rownames(HR.matrix) = HR.matrix[,1]
HR.matrix = HR.matrix[,-c(1)]
mode(HR.matrix) = "numeric"

HR.table$p_value = signif(HR.table$p_value, 2)
HR.table$HR = signif(HR.table$HR, 3)
tabletext<-cbind(
  c("IMS", as.character(HR.table$IMS)[c(1:n.cells)]),
  c("p-value", HR.table$p_value[c(1:n.cells)]),
  c("HR",      HR.table$HR[c(1:n.cells)]),
  c("N", HR.table$N[c(1:n.cells)]))


dir.create("./Figures/C39.4_Forest_plot", showWarnings = FALSE)
pdf(file = paste0("./Figures/C39.4_Forest_plot/C39.4_Forest_plot",
                  "_", Ethnicity,".pdf"),
    height = 5, width = 4)

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
           boxsize = .10,
           vertices = FALSE,
           col=fpColors(box="darkblue",line="darkgrey"),
           txt_gp = fpTxtGp(label = gpar(fontsize = 7), xlab = gpar(fontsize = 7),
                            ticks = gpar(fontsize = 10))
)
dev.off()


