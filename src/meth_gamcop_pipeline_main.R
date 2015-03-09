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
source("src/visualize_diff_exp_dlg.R")
source("src/meth_functions.R")
source("src/gene_list_overlap.R")

######################################################################################################################
# reading in the data
setwd("/Volumes/StorageDisk/Meth_DF5/pipeline")
# Clinical data file
#clinMat <- readFeatureMatrix("basis/data_CLIN_20150203.fm")
clinMat <- readFeatureMatrix("basis/data_CLIN_20150303.fm")
#methMat <- readFeatureMatrix("products/data_METH_20140129_norm_filtered_outlier_logit_admix_limma.fm")
load(file="/Volumes/StorageDisk/Meth_DF5/pipeline/methmat_feb_9.rda")
#batchMat <- read.csv("/Volumes/StorageDisk/Meth_DF5/pipeline/basis/DF5_Methylation_Batches.csv", stringsAsFactors=FALSE)
#clinMat <- methAddBatch(clinMat, batchMat, "N:M:METH:Data:MethBatch")

idx <- grep(rownames(clinMat), pattern = "Preec")
rownames(clinMat)[idx] <- "B:M:CLIN:Critical_Phenotype:Preeclampsia_Eclampsia"

bloodVar <- "N:M:CLIN:Data:Date_of_Blood_Collection_20150223__relative_to_Date_of_Birth"
#bloodVar <- "N:M:CLIN:Data:Date_of_Blood_Collection__relative_to_Date_of_Birth"  

bloodDrawDates <- as.numeric(clinMat[bloodVar, ])

targets <- rownames(clinMat)[grep(pattern="Critical_Phenotype", rownames(clinMat))]

dirs <- c("/Volumes/StorageDisk/Meth_DF5/pipeline/boot_DE_blood_all_days/",
          "/Volumes/StorageDisk/Meth_DF5/pipeline/boot_DE_blood_day_04/",
          "/Volumes/StorageDisk/Meth_DF5/pipeline/boot_DE_blood_day_1/")

dayString <- c("AllDays", "04Days", "1Days")


for (day in 1:3) {
  if (day == 1) {
    # then take them all 
    idx <- bloodDrawDates > -100000
    covariates <- c("N:M:SURV:Data:Date_of_Birth__relative_to_Date_of_Birth", 
                    "C:M:ADMX:Data:Admix_80_Percent",
                    bloodVar)
  } else if (day == 2) {
    idx <- bloodDrawDates >= 0 & bloodDrawDates <= 4
    covariates <- c("N:M:SURV:Data:Date_of_Birth__relative_to_Date_of_Birth", 
                    "C:M:ADMX:Data:Admix_80_Percent",
                    bloodVar)
  } else if (day == 3) {
    idx <- bloodDrawDates == 1
    covariates <- c("N:M:SURV:Data:Date_of_Birth__relative_to_Date_of_Birth", 
                    "C:M:ADMX:Data:Admix_80_Percent")
  }
  idx[is.na(idx)] <- F
  clinMatFilt <- clinMat[,idx]  
  outdir <- dirs[day]
  
    for (ta in targets) {
      try({
        targetString <- str_split(ta, ":")[[1]][5]; print(targetString)
        #Rprof("boot_profile.out")
        Rprof("boot_reg_profile.out")
        system.time(deTable <- bootDiffFun(clinMat=clinMatFilt, dataMat=methMat, targetPheno=ta,
                                covarVec=covariates, FCThresh=0.001, pValueThresh=0.05,
                                writingDir=outdir, reps=10, cpus=4, robustFlag=F))
        Rprof(NULL)
        summaryRprof("boot_reg_profile.out")
        proftable("boot_reg_profile.out")
        #deTable <- diffExprFun(clinMat=clinMatFilt, dataMat=methMat, targetPheno=ta,
        #                       covarVec=covariates, FCThresh=0.001, pValueThresh=0.05,
        #                       writingDir=outdir, robustFlag=T)
        print("*****************")
        print(day)
        print(targetString)
        print(sum(clinMatFilt[ta,] == "true", na.rm = T))
        print(sum(clinMatFilt[ta,] == "false", na.rm = T))
        print(ncol(clinMatFilt[ta,]))
        print(sum(clinMatFilt[ta,] == "true", na.rm = T)/ncol(clinMat[ta,]))
        print(sum(deTable$adj.P.Val < 0.05, na.rm = T))
        print("*****************")
      
        geneTable <- mapToGenes(deTable = deTable, 1.010)
        ddd <- cbind(deTable, geneTable)
        foutstring <- add_date_tag(str_join(outdir, "METH_DE_", dayString[day], "_", targetString), ".txt")
        write.table(ddd, quote=F, file=foutstring, sep="\t")
      }) 
    }
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

