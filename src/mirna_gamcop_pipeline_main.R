# Script to run miRNA-Seq differential expression analysis using the unified GAMCOP Pipeline
# Date: 02/02/2015
# Author: Nyasha Chambwe
# Source functions from: https://github.com/Gibbsdavidl/gamcop_pipeline.git
# Developed in collaboration between Crystal Humphries, Dave Gibbs and Nyasha Chambwe
source("src/diffExprFun.R")
source("src/diffExpMeth.R")
source("src/readFeatureMatrix.R")
source("src/functions_visualize_diff_exp.R")
source("src/visualize_diff_exp.R")

# Clinical data file
clinDat <-  read.table("/Volumes/~nchambwe-1/inova/projects/gamcop/smrna-analysis/data/clinical/data_CLIN_20150128.fm", check.names=F, header=T, row.names = 1)
 
molDat <- read.table("/Volumes/~nchambwe-1/inova/projects/gamcop/smrna-analysis/data/molecular/small_rna_seq/Study_101/batch_summaries/data_MIRSeq_LIMMA_corrected_20141201.fm", 
                                header=T, sep="\t", row.names=1, check.names = F)
#covariates <- c("N:M:SURV:Data:Date_of_Birth__relative_to_Date_of_Birth", "N:M:SURV:Data:Years_Drinking")
covariates <- c("B:M:CLIN:Data:Vaginal_Bleeding_Spotting", "B:M:CLIN:Data:Prom")
#targetPhenotype <- "B:NB:CLIN:Critical_Phenotype:Preterm"
targetPhenotype <- "N:NB:CLIN:Critical_Phenotype:TermCategory"
#B:NB:CLIN:Critical_Phenotype:Preterm
#N:NB:CLIN:Critical_Phenotype:TermCategory
#N:NB:CLIN:Critical_Phenotype:Gestational_Age_at_Delivery
#B:NB:MRGE:Critical_Phenotype:Immune_Related
diffExprFun(clinMat = clinDat, dataMat=molDat, targetPheno=targetPhenotype, covarVec=covariates,FCThresh=0.05, pValueThresh=0.1, 
            writingDir="/Volumes/~nchambwe-1/inova/projects/gamcop/") 
    