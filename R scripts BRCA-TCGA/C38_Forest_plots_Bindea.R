
# Setup environment
rm(list = ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA")  

source("../tools/ipak.function.R")
required.packages <- c("RColorBrewer", "forestplot")
ipak(required.packages)

# Set parameters
Outcome = "OS"
Ethnicity = "All"  # "All", "White", "Black"
IMS = "All" # "BasalMyo"

# Load
load(paste0("./Analysis/C37_Coxph_Bindea/Coxph_table_", IMS, "_", Ethnicity, ".Rdata"))

HR.table = results[-which(results$Cell_types == "Tgd"),]
HR.table = HR.table[-which(HR.table$Cell_types == "Mast cells"),]

n.cells = nrow(HR.table)
x = n.cells + 2

if(Ethnicity == "All" & IMS == "All"){
  HR.table = HR.table[order(HR.table$HR),]
  order_cells = as.character(HR.table$Cell_types)
  save(order_cells, file = "./Analysis/C37_Coxph_Bindea/Order_cells.Rdata")
}
if(Ethnicity == "All" & IMS == "BasalMyo"){
  HR.table = HR.table[order(HR.table$HR),]
  order_cells = as.character(HR.table$Cell_types)
  save(order_cells, file = "./Analysis/C37_Coxph_Bindea/Order_cells_BasalMyo.Rdata")
}
if(Ethnicity %in% c("Black", "White") & IMS == "All"){
  load("./Analysis/C37_Coxph_Bindea/Order_cells.Rdata")
  HR.table = HR.table[order(match(HR.table$Cell_types, order_cells)),]
} 
if(Ethnicity %in% c("Black", "White") & IMS == "BasalMyo"){
  load("./Analysis/C37_Coxph_Bindea/Order_cells_BasalMyo.Rdata")
  HR.table = HR.table[order(match(HR.table$Cell_types, order_cells)),]
}


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


dir.create("./Figures/C38_Forest_plot_Bindea", showWarnings = FALSE)
pdf(file = paste0("./Figures/C38_Forest_plot_Bindea/C38_new_Forest_plot_Bindea",
                  IMS, "_", Ethnicity,".pdf"),
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
           boxsize = .25,
           vertices = FALSE,
           col=fpColors(box="darkblue",line="darkgrey"),
           txt_gp = fpTxtGp(label = gpar(fontsize = 7), xlab = gpar(fontsize = 7),
                            ticks = gpar(fontsize = 10))
           )
dev.off()


