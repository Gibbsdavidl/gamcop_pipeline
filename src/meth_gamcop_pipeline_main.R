# Script to run METH differential expression analysis using the unified GAMCOP Pipeline
# Date: 02/04/2015
# Author: Nyasha Chambwe + Modified by David Gibbs
# Source functions from: https://github.com/Gibbsdavidl/gamcop_pipeline.git
# Developed in collaboration between Crystal Humphries, Dave Gibbs and Nyasha Chambwe

setwd("/Volumes/StorageDisk/Meth_DF5/gamcop_pipeline")
source("src/differential.R")
source("src/diffExprFun.R")
source("src/readFeatureMatrix.R")
source("src/functions_visualize_diff_exp.R")
source("src/visualize_diff_exp.R")
source("src/meth_functions.R")
source("src/gene_list_overlap.R")

######################################################################################################################
# reading in the data
setwd("/Volumes/StorageDisk/Meth_DF5/pipeline")
# Clinical data file
clinMat <- readFeatureMatrix("basis/data_CLIN_20150203.fm")
#methMat <- readFeatureMatrix("products/data_METH_20140129_norm_filtered_outlier_logit_admix_limma.fm")
load(file="/Volumes/StorageDisk/Meth_DF5/pipeline/methmat_feb_9.rda")
#batchMat <- read.csv("/Volumes/StorageDisk/Meth_DF5/pipeline/basis/DF5_Methylation_Batches.csv", stringsAsFactors=FALSE)
#clinMat <- methAddBatch(clinMat, batchMat, "N:M:METH:Data:MethBatch")

#####################################################################################################################
# fixing clinmat
idx <- grep(rownames(clinMat), pattern = "Preec")
rownames(clinMat)[idx] <- "B:M:CLIN:Critical_Phenotype:Preeclampsia_Eclampsia"
bloodDrawDates <- as.numeric(clinMat["N:M:CLIN:Data:Date_of_Blood_Collection__relative_to_Date_of_Birth", ])
#idx <- bloodDrawDates == 1
#idx <- bloodDrawDates >= 0 & bloodDrawDates <= 4
idx[is.na(idx)] <- F
clinMat <- clinMat[,idx]
#######################################################################################################################

covariates <- c("N:M:SURV:Data:Date_of_Birth__relative_to_Date_of_Birth", 
                "C:M:ADMX:Data:Admix_80_Percent")

covariates <- c("N:M:SURV:Data:Date_of_Birth__relative_to_Date_of_Birth", 
                "C:M:ADMX:Data:Admix_80_Percent",
                "N:M:CLIN:Data:Date_of_Blood_Collection__relative_to_Date_of_Birth")

targets <- c(
"B:NB:CLIN:Critical_Phenotype:Preterm",
"N:NB:CLIN:Critical_Phenotype:TermCategory",
"B:NB:CLIN:Critical_Phenotype:1n2v4",
"B:NB:CLIN:Critical_Phenotype:1v4",
"B:NB:MRGE:Critical_Phenotype:History_of_Preterm_Birth",
"B:NB:MRGE:Critical_Phenotype:Hypertension_Related",
"B:NB:MRGE:Critical_Phenotype:Immune_Related",
"B:M:CLIN:Critical_Phenotype:Incompetent_Shortened_Cervix",
"B:M:CLIN:Critical_Phenotype:Inova_Idiopathic_NA",
"B:M:CLIN:Critical_Phenotype:Preeclampsia_Eclampsia",
"B:NB:MRGE:Critical_Phenotype:Uterine_Related")

for (ta in targets) {
  targetString <- str_split(ta, ":")[[1]][5]; print(targetString)
  deTable <- diffExprFun(clinMat=clinMat, dataMat=methMat, targetPheno=ta,
                covarVec=covariates, FCThresh=0.001, pValueThresh=0.05,
                writingDir="/Volumes/StorageDisk/Meth_DF5/pipeline/DE_blood_all_days") 
  geneTable <- mapToGenes(deTable = deTable, 1.010)
  ddd <- cbind(deTable, geneTable)
  foutstring <- str_join("/Volumes/StorageDisk/Meth_DF5/pipeline/DE_blood_all_days/", targetString, ".txt", collapse = "")
  write.table(ddd, quote=F, file=foutstring, sep="\t")
}

#############################################################################################################################
#Overlap tests
dbPTB_Genes_Feb11_2013 <- read.csv("/Volumes/StorageDisk/Meth_DF5/pipeline/basis/dbPTB_Genes_Feb11_2013.csv", header=FALSE, stringsAsFactors=FALSE)
universe <- makeUniverse()
dbPTB <- dbPTB_Genes_Feb11_2013$V1[dbPTB_Genes_Feb11_2013$V1 %in% universe]

outputFiles <- c("1n2v4.txt", "1v4.txt", "Preterm.txt",
  "History_of_Preterm_Birth.txt", "Preeclampsia_Eclampsia.txt",
  "Hypertension_Related.txt", "Immune_Related.txt",
  "Incompetent_Shortened_Cervix.txt", "Inova_Idiopathic_NA.txt",
  "Preeclampsia_Eclampsia.txt", "TermCategory.txt", "Uterine_Related.txt")

for (outf in outputFiles) {
  fff <- str_join("/Volumes/StorageDisk/Meth_DF5/pipeline/DE_blood_day1/",outf,collapse="")
  dat <- read.delim(fff, stringsAsFactors=FALSE)
  diffExprGenes <- unique(as.character(dat$nearestGeneSymbol[dat$adj.P.Val <= 0.05]))
  if (length(diffExprGenes) > 0) {
    print(outf)
    print(hypergeomTest(diffExprGenes, dbPTB, universe))
  }
}

#############################################################################################################################

