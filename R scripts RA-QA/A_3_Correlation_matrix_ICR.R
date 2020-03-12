#######
#
# Correlation matrix between genes of interest
#
#######


# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-RA-QA")
source("../tools/ipak.function.R")
required.packages <- c("corrplot", "stringr")
ipak(required.packages)  

# Set parameters
selected_genes = "ICR"
test = "pearson"
colpattern = colorRampPalette(c("#152B7E", "white", "#1B7E09"))(n = 297)

# load data
load("./Data/expression matrix/CTA_TNBC_BRCA_Qatar_RNASeq.tumor.tissue.filtered.EDAseq.QN.Rdata")

if (selected_genes == "ICR") {
  load("../tools/ICR_genes.RData")
  genes_to_correlate = ICR_genes
}

subset_RNAseq_log2 = t(RNASeq.QN.LOG2[row.names(RNASeq.QN.LOG2) %in% genes_to_correlate, ])

# Correlation matrix
ICR_cor <- cor (subset_RNAseq_log2,method=test)

# Correlation significance
cor.mtest <- function(mat, conf.level = 0.95) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], method = test, conf.level = conf.level)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      lowCI.mat[i, j] <- lowCI.mat[j, i] <- tmp$conf.int[1]
      uppCI.mat[i, j] <- uppCI.mat[j, i] <- tmp$conf.int[2]
    }
  }
  return(list(p.mat, lowCI.mat, uppCI.mat))
}
ICR_cor_sign <- cor.mtest(subset_RNAseq_log2, 0.95)

# Correlation plot
dir.create("./Figures/Correlation_plots", showWarnings = FALSE)
png(paste0("./Figures/Correlation_plots/", selected_genes, "_Correlation_plot.png"),
    res=600,height=6,width=6,unit="in")

cex.before <- par("cex")
par(cex = 0.45)
lims=c(-1,1)
if (length(ICR_cor[ICR_cor<0]) == 0) {lims=c(0,1)}
annotation = data.frame (gene = rownames(ICR_cor),color = c(rep("#CC0506",20)),stringsAsFactors = FALSE)
annotation$color[annotation$gene %in% c("IDO1","CD274","CTLA4","FOXP3","PDCD1")] = "#41719C"
annotation = annotation[corrMatOrder(ICR_cor,order="FPC"),]

mean_correlation = round(mean(ICR_cor),2)
corrplot.mixed (ICR_cor,
                #type="lower",
                #p.mat = ICR_cor_sign[[1]],                                                                      # add significance to correlations
                lower.col = colpattern,
                upper.col = colpattern,
                lower = "square",
                upper ="number",
                order="FPC",
                cl.lim=lims,                                                                                               # only positive correlations
                tl.pos ="lt",
                tl.col = as.character(annotation$color),
                insig= "pch",                                                                                              # remove insignificant correlations
                pch = "x",
                pch.cex= 1.5,
                tl.cex = 1/par("cex"),
                cl.cex = 1/par("cex"),
                cex.main = 1/par("cex"),
                mar=c(6,4.1,7,5))
title(main = list(paste0("Correlation between ", selected_genes, " genes. \n ","Mean: ", mean_correlation,". Number of patients: ", nrow(subset_RNAseq_log2), "."),
                  cex = 2.2), line = -2.5)
title(sub = list(paste0("Significance level of correlation is represented by the size of the squares."), cex = 1.5), line = 0)
par(cex = cex.before)
dev.off()
