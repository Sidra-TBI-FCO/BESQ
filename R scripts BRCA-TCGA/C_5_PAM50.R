#######
#
# Correlation matrix between genes of interest
# PAM50 obtained from: https://www.biostars.org/p/210106/
#
#######


# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/NNN-BRCA/RNAseq-TCGA_BRCA")
source("../tools/ipak.function.R")
required.packages <- c("heatmap.plus", "impute", "ctc", "amap")
ipak(required.packages)  

#source("./Analysis QBRI and TCGA/R tools/PAM50_Parker/")

#Set parameters

# Load data
load("../tools/ICR_genes.RData")
load("./Analysis/ICR data/TCGA_BRCA_table_cluster_assignment.RData")
Clinical_data = read.csv("./Data/Clinical Data/TCGA_CLINICAL_DATA_CELL_2018_S1.csv", stringsAsFactors = FALSE)
load("./Data/expression matrix/BRCA_gene_RNAseq_normalized_Tissue_Type_and_whitelists_filtered.Rdata")

# PAM50 classification
## RNASeq Data from the EDASeq protocol after log2 transformation
RNASeq.QN.LOG2 = log(filtered.norm.RNAseqData +1, 2)

#rename to old genenames to fit with this method
rownames(RNASeq.QN.LOG2)[rownames(RNASeq.QN.LOG2) == "NDC80"] = "KNTC2"
rownames(RNASeq.QN.LOG2)[rownames(RNASeq.QN.LOG2) == "NUF2"] = "CDCA1"
rownames(RNASeq.QN.LOG2)[rownames(RNASeq.QN.LOG2) == "ORC6"] = "ORC6L"

dir.create("./Analysis/PAM50", showWarnings = FALSE)

write.table(RNASeq.QN.LOG2,file= "./Analysis/PAM50/RNASeq_QN_LOG2_Input_for_PAM50.txt",sep = "\t",quote=FALSE,col.names=NA)

paramDir<- "../tools/PAM50_Parker/bioclassifier_R"
inputDir<- "./Analysis/PAM50"        # the location of the data matrix, and where output will be located

inputFile<- "RNASeq_QN_LOG2_Input_for_PAM50.txt" # the input data matrix as a tab delimited text file
short<-"TCGA_BRCA_IMS_Output"       # short name that will be used for output files

calibrationParameters<- NA 	              #the column of the "mediansPerDataset.txt" file to use for calibration; 
#NA will force centering within the test set & -1 will not do any 
#adjustment (when adjustment performed by used)

hasClinical<-FALSE 	#may include tumor size as second row, with 'T' as the gene name, 
#and encoded as binary (0 for size <= 2cm or 1 for size > 2cm)
#set this variable to FALSE if tumor size is not available

collapseMethod<-"mean" # can be mean or iqr (probe with max iqr is selected)
# typically, mean is preferred for long oligo and
# iqr is preferred for short oligo platforms


####
# run the assignment algorithm
####

y = unlist(y)

source(paste(paramDir,"subtypePrediction_functions.R",sep="/"))
source(paste(paramDir,"subtypePrediction_distributed.R",sep="/"))

y <- data.frame(matrix(unlist(y), nrow=6, byrow=T))


#filterout PAM50 to samples needed 
TCGA.IMS <- read.table ("./Analysis/PAM50/TCGA_BRCA_IMS_Output_pam50scores.txt",header=TRUE)

#row.names.remove <- c("05.06", "09.05", "10.24") # In Aisha's thesis removed because of other ethnicity: Iran, Canada, Irish

rownames(TCGA.IMS) <- TCGA.IMS[,1]
TCGA.IMS <- TCGA.IMS[,-1]

# Add information to clinical data
Clinical_data = Clinical_data[which(Clinical_data$type == "BRCA"),]
Clinical_data$IMS = TCGA.IMS$Call[match(Clinical_data$bcr_patient_barcode, substring(rownames(TCGA.IMS), 1, 12))]

Clinical_data_ann = Clinical_data
save(Clinical_data_ann, file = "./Analysis/Sample_annotations.Rdata")
