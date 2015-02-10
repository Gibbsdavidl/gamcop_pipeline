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

setwd("/Volumes/StorageDisk/Meth_DF5/pipeline")

# Clinical data file
clinMat <- readFeatureMatrix("basis/data_CLIN_20150203.fm")
methMat <- readFeatureMatrix("products/data_METH_20140129_norm_filtered_outlier_logit_admix_limma.fm")
#batchMat <- read.csv("/Volumes/StorageDisk/Meth_DF5/pipeline/basis/DF5_Methylation_Batches.csv", stringsAsFactors=FALSE)
#clinMat <- methAddBatch(clinMat, batchMat, "N:M:METH:Data:MethBatch")


bloodDrawDates <- as.numeric(clinMat["N:M:CLIN:Data:Date_of_Blood_Collection__relative_to_Date_of_Birth", ])
idx <- bloodDrawDates >= 0 & bloodDrawDates <= 1
idx[is.na(idx)] <- F
clinMat <- clinMat[,idx]


covariates <- c("N:M:SURV:Data:Date_of_Birth__relative_to_Date_of_Birth", "C:M:ADMX:Data:Admix_80_Percent")
targetPhenotype <- "B:NB:CLIN:Critical_Phenotype:Preterm"
#targetPhenotype <- "N:NB:CLIN:Critical_Phenotype:TermCategory"
#"N:CLIN:Date_of_Blood_Collection__relative_to_Date_of_Birth:M::::"
#B:NB:CLIN:Critical_Phenotype:Preterm
#N:NB:CLIN:Critical_Phenotype:TermCategory
#N:NB:CLIN:Critical_Phenotype:Gestational_Age_at_Delivery

deTable <- diffExprFun(clinMat=clinMat, dataMat=methMat, targetPheno=targetPhenotype,
              covarVec=covariates, FCThresh=0.001, pValueThresh=0.05,
              writingDir="/Volumes/StorageDisk/Meth_DF5/pipeline/DE") 

geneTable <- mapToGenes(deTable = deTable, 0.05)
write.table(geneTable, quote=F, file="/Volumes/StorageDisk/Meth_DF5/pipeline/DE/termcat_genes_mothers_0-1.txt", sep="\t")


##############################################

#Overlap tests
dbPTB_Genes_Feb11_2013 <- read.csv("/Volumes/StorageDisk/Meth_DF5/pipeline/basis/dbPTB_Genes_Feb11_2013.csv", header=FALSE, stringsAsFactors=FALSE)

universe <- read.table("/Volumes/StorageDisk/Meth_DF5/pipeline/products/Filtered_Meth_Probe_Gene_Universe.txt", quote="\"", stringsAsFactors=FALSE)

unique(pretermGeneTable$nearestGeneSymbol)
universe <- makeUniverse()
dbPTB <- dbPTB_Genes_Feb11_2013$V1[dbPTB_Genes_Feb11_2013$V1 %in% universe]
deGenes <- unique(as.character(termcat04GeneTable$nearestGeneSymbol))

a <- length(intersect(deGenes, dbPTB))  # sig and in dbPTB
b <- length(deGenes) - a   # sig and not in dbPTB
c <- length(dbPTB) - a # non sig and IN dbPTB
d <- length(intersect(setdiff(universe, deGenes), setdiff(universe, dbPTB)))
a
b
c
d
a+b
c+d

phyper(a, length(dbPTB), (length(universe) - length(dbPTB)), length(deGenes), lower.tail = T)
# Numerical parameters in order:
# (success-in-sample, success-in-bkgd, failure-in-bkgd, sample-size).

#phyper(x,m,n,k)
#x, vector of quantiles representing the number of white balls drawn without replacement from an urn which contains both black and white balls.
#m	 the number of white balls in the urn.
#n	 the number of black balls in the urn.
#k	 the number of balls drawn from the urn.





