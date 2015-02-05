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


covariates <- c("N:M:SURV:Data:Date_of_Birth__relative_to_Date_of_Birth", "C:M:ADMX:Data:Admix_80_Percent")
targetPhenotype <- "B:NB:CLIN:Critical_Phenotype:Preterm"
#"N:CLIN:Date_of_Blood_Collection__relative_to_Date_of_Birth:M::::"
#targetPhenotype <- "N:NB:CLIN:Critical_Phenotype:TermCategory"
#B:NB:CLIN:Critical_Phenotype:Preterm
#N:NB:CLIN:Critical_Phenotype:TermCategory
#N:NB:CLIN:Critical_Phenotype:Gestational_Age_at_Delivery

diffExprFun(clinMat=clinMat, dataMat=methMat, targetPheno=targetPhenotype,
            covarVec=covariates, FCThresh=0.263, pValueThresh=0.05,
            writingDir="/Volumes/StorageDisk/Meth_DF5/pipeline/DE") 

