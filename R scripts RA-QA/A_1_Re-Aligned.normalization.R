#script for processing and normalistaton of RNAseq data from subreads
# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/External Collaborations/NNN-BRCA/RNAseq Data/expression matrix")
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-RA-QA/Data/expression matrix")
source("~/Dropbox (TBI-Lab)/DB-LAB/bioinformatics tools/R scripts/ipak.function.R")

#read text file
input.file = "./counts.txt"
subread.command <- readLines (input.file,n = 1)
raw.data = read.csv (input.file,stringsAsFactors = FALSE,skip = 1,sep = "\t")

#extract relevant data
gene.data <- raw.data[,c(1:6)]
expression.data <- raw.data[,-c(2:6)]
#expression.data <- expression.data[,-ncol(expression.data)]
rownames(expression.data) <- expression.data$Geneid
expression.data$Geneid <- NULL
colnames(expression.data) <- gsub(x = colnames(expression.data),pattern = ".bam",replacement = "")
colnames(expression.data) <- gsub(x = colnames(expression.data),pattern = "BAM.",replacement = "")
colnames(expression.data) <- gsub(x = colnames(expression.data),pattern = "CTA_TNBC_BRCA_Qatar_",replacement = "")
colnames(expression.data) <- gsub(x = colnames(expression.data),pattern = ".trimmed",replacement = "")
expression.matrix <- as.matrix(expression.data)
#save as R datafile
save (gene.data,subread.command,expression.matrix,file="CTA_TNBC_BRCA_Qatar_RNASeq.dataset.HPC.Rdata")

# normalizationa
# dependencies
ibiopak("EDASeq")
ibiopak("base64enc")
ibiopak("preprocessCore")
load ("~/Dropbox (TBI-Lab)/DB-LAB/bioinformatics tools/GeneInfo/geneInfo.July2017.RData")                                                               

# cleanup expression data
geneInfo <- as.data.frame(geneInfo)
geneInfo <- geneInfo[!is.na(geneInfo[,1]),]                                                                          # drop the genes without information and convert to data frame
available.genes <- unique(rownames(expression.matrix)[which(rownames(expression.matrix) %in% rownames(geneInfo))])   # genes in geneinfo and RNASeqdata
geneInfo <- geneInfo[available.genes,]                                                                               # drop the GCcontent info on unavailable genes
expression.filtered <- expression.matrix[available.genes,]                                                           # drop the genes without information from the RNAseq.DATA
mode(expression.filtered) <- "numeric"
dim(expression.filtered)

save(expression.filtered, file = "./CTA_TNBC_BRCA_Qatar_RNASeq.filtered.Rdata")

#split tumor tissue and cell lines
expression.filtered.tissue <- expression.filtered[,c(1:28)]
expression.filtered.cells <- expression.filtered[,-c(1:28)]
colnames(expression.filtered.tissue) <- paste0("S.",colnames(expression.filtered.tissue))

#annotation
cell.line.samples <- data.frame(sample_ID=colnames(expression.filtered.cells),
                                cell_line=c("468","468","BT549","BT549","BT549","HCC1954","HCC1954","HCC1954","HCC70","HCC70","MDA.MB.231","MDA.MB.231"),
                                siRNA=c("NC","PRAME","siLDHC","siPRAME","siCTRL","siCTRL","siLDHC","siPRAME","siCTRL","siPRAME","siCTRL","siLDHC")
)

# QC
dim (expression.filtered.tissue)
dim (expression.filtered.cells)
#remove flat genes
#expression.filtered.tissue <- expression.filtered.tissue[-which(rowMeans(expression.filtered.tissue)==0),] #708
#dim(expression.filtered.tissue)
#expression.filtered.cells <- expression.filtered.cells[-which(rowMeans(expression.filtered.cells)==0),]    #1230
#dim(expression.filtered.cells)
#remove not finite values
all(is.finite(expression.filtered.tissue)) #TRUE
all(is.finite(expression.filtered.cells))  #TRUE

#total count stats
dev.new
par(mar=c(15,4,4,2))
barplot(colSums(expression.filtered.cells),las=2)
barplot(colSums(expression.filtered.tissue))

#PCA
expression.filtered.pca <- expression.filtered.cells
pca <- prcomp(t(expression.filtered.pca), center= TRUE, scale=TRUE)
scores <- data.frame(pca$x[,1:3])
scores <- scores[cell.line.samples$sample_ID,]
scores$label <- cell.line.samples$cell_line
scores$shape <- as.character(cell.line.samples$siRNA)
scores$shape[scores$shape=="siCTRL"] <- 16
scores$shape[scores$shape=="siLDHC"] <- 17
scores$shape[scores$shape=="siPRAME"] <- 18
scores$shape[scores$shape=="NC"] <- 15
scores$shape[scores$shape=="PRAME"] <- 14
scores$shape <- as.numeric(scores$shape)
myColors <- data.frame(color=brewer.pal(length(unique(scores$label)),"Set1"))
rownames(myColors) <- levels (as.factor(scores$label))
scores$color  <- myColors$color[match(scores$label,rownames(myColors))]
#3Dplot
dev.new()
s3d <- scatterplot3d(scores[,c(1:3)],
                     #pch = 19,
                     type="h",
                     cex.symbols=1.5,
                     color=scores$color,
                     pch = scores$shape,
                     main="3-D PCA cell lines",
                     xlab="PC1 (X%)",
                     ylab="PC2 (Y%)",
                     zlab="PC3 (Z%)"
)
text(s3d$xyz.convert(scores[,c(1:3)]), labels = rownames(scores),
     cex= 2, col = "steelblue")
legend("bottom", legend = c(rownames(myColors),"",levels(cell.line.samples$siRNA)),
       col = c(as.character(myColors$color),"white",rep("grey",5)),
       pch = c(16,16,16,16,16,16,15,14,16,17,18),
       inset = -0.15,
       xpd = TRUE,
       horiz = TRUE,
       #lty= 1,
       #lwd = 5,
       cex = 1.5)


# EDAseq cell lines
RNASeq.expr.set <- newSeqExpressionSet(expression.filtered.cells , featureData = geneInfo)                     # create a new SeqExpressionSet object.
fData(RNASeq.expr.set)[, "gcContent"] <- as.numeric(geneInfo[, "gcContent"])                             # make sure gcContenet is numeric
RNASeq.expr.set <- withinLaneNormalization(RNASeq.expr.set, "gcContent", which = "upper", offset = TRUE) # removes lane gene specific effects, for example effects related to gene length or GC content
RNASeq.expr.set <- betweenLaneNormalization(RNASeq.expr.set, which = "upper", offset = TRUE)             # removes effect related to in between lane distributional differences, as sequencing depth
RNASeq.NORM <-  log(expression.filtered.cells + .1) + offst(RNASeq.expr.set)                                   # apply the Edaseq Ofset
RNASeq.NORM <-  floor(exp(RNASeq.NORM) - .1)                                                             # return non decimal values

# Quantile normalisation RNA
RNASeq.NORM.quantiles <- normalize.quantiles(RNASeq.NORM)                                                # Quantile normailze
RNASeq.NORM.quantiles <- floor(RNASeq.NORM.quantiles)                                                    # return non decimal values
rownames(RNASeq.NORM.quantiles) <- rownames(RNASeq.NORM)
colnames(RNASeq.NORM.quantiles) <- colnames(RNASeq.NORM)
# log transformation
RNASeq.QN.LOG2<-log(RNASeq.NORM.quantiles+1,2)
save (RNASeq.QN.LOG2,RNASeq.NORM.quantiles,geneInfo,
      file="./CTA_TNBC_BRCA_Qatar_RNASeq.cell.lines.filtered.EDAseq.QN.Rdata")
#plot Log2 data
dev.new()
ggplot(melt(RNASeq.QN.LOG2)) +
  geom_freqpoly(aes(x = value,y = ..density.., colour = Var2)) +
  xlim(c(2,20))


# EDAseq tissue
expression.filtered.tissue <- expression.filtered.tissue[,colSums(expression.filtered.tissue) > 50000]
save(expression.filtered.tissue, file = "./RNASeq.filtered.tissue.raw.Rdata")
RNASeq.expr.set <- newSeqExpressionSet(expression.filtered.tissue , featureData = geneInfo)              # create a new SeqExpressionSet object.
fData(RNASeq.expr.set)[, "gcContent"] <- as.numeric(geneInfo[, "gcContent"])                             # make sure gcContenet is numeric
RNASeq.expr.set <- withinLaneNormalization(RNASeq.expr.set, "gcContent", which = "upper", offset = TRUE) # removes lane gene specific effects, for example effects related to gene length or GC content
RNASeq.expr.set <- betweenLaneNormalization(RNASeq.expr.set, which = "upper", offset = TRUE)            # removes effect related to in between lane distributional differences, as sequencing depth
RNASeq.NORM <-  log(expression.filtered.tissue + .1) + offst(RNASeq.expr.set)                            # apply the Edaseq Ofset
RNASeq.NORM <-  floor(exp(RNASeq.NORM) - .1)                                                             # return non decimal values
m<-counts(RNASeq.expr.set)
all(is.finite(m))

# Quantile normalisation RNA
RNASeq.NORM.quantiles <- normalize.quantiles(RNASeq.NORM)                                                # Quantile normailze
RNASeq.NORM.quantiles <- floor(RNASeq.NORM.quantiles)                                                    # return non decimal values
rownames(RNASeq.NORM.quantiles) <- rownames(RNASeq.NORM)
colnames(RNASeq.NORM.quantiles) <- colnames(RNASeq.NORM)
# log transformation
RNASeq.QN.LOG2<-log(RNASeq.NORM.quantiles+1,2)
save (RNASeq.QN.LOG2,RNASeq.NORM.quantiles,geneInfo,
      file="./CTA_TNBC_BRCA_Qatar_RNASeq.tumor.tissue.filtered.EDAseq.QN.Rdata")
#plot Log2 data
dev.new()
ggplot(melt(RNASeq.QN.LOG2)) +
  geom_freqpoly(aes(x = value,y = ..density.., colour = Var2)) +
  xlim(c(2,20))
dev.new()
expression.filtered.tissue
ggplot(melt(log(expression.filtered.tissue+0.1))) +
  geom_freqpoly(aes(x = value,y = ..density.., colour = Var2)) +
  xlim(c(2,20))
