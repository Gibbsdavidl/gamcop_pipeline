# Script to run mRNA-Seq differential expression analysis using the unified GAMCOP Pipeline
# Date: 03/09/2015
# Author: Crystal Humphries  (derived from Dave's Methylation Pipeline)

######################################################################################################################

# Source Functions and R Packages -----------------------------------------
setwd("/Users/chumphri/Projects/RNAseq/Analysis_in_R/gamcop_pipeline")

# Source functions from: https://github.com/Gibbsdavidl/gamcop_pipeline.git
# Developed in collaboration between Crystal Humphries, Dave Gibbs and Nyasha Chambwe
source("src/diffExprFun.R")
source("src/differential.R")
source("src/readFeatureMatrix.R")
source("src/functions_visualize_diff_exp.R")
source("src/visualize_diff_exp.R")
source("src/additional_data_parsing.R")

# Returns limma diffExp table of signficant hits according to pval and fc thresholds
# Using names may contribute to misaligned headers
filter_diffExp <- function(pathToDeTable, FCThresh, pValueThresh){
  tab <- read.table(pathToDeTable, header = T, sep="")
  stopifnot( "logFC" %in% names(tab), "adj.P.Val" %in% names(tab))
  diffExpTab <- tab[tab$adj.P.Val <= pValueThresh & abs(tab$logFC) >= FCThresh,]
  return(diffExpTab)
}

######################################################################################################################

# Load Input Files --------------------------------------------------------
# Load data
setwd("/Users/chumphri/Projects/RNAseq/Feature Matrices/")

# Clinical data file
clinDat <-  read.table("data_CLIN_20150303.fm", check.names=F, header=T, row.names = 1, stringsAsFactors = F)

molDat <- read.table("data_RNASeq_outliers_removed_normalized_GenesWithCountsAbove2_batchCorrected.fm", header=T, sep="\t", row.names=1, check.names = F)

colnames(clinDat)<-cor_colnames(clinDat)
colnames(molDat) <- cor_colnames(molDat)

# subset clinical matrix for samples with available molecular data
clinDat <- clinDat[,(colnames(clinDat)  %in% colnames(molDat))]

# Fix forward slash in Preeclampsia/Eclampsia
idx <- grep(rownames(clinDat), pattern = "Preec")
rownames(clinDat)[idx] <- "B:M:CLIN:Critical_Phenotype:Preeclampsia_Eclampsia"

targetPhenotypes <- row.names(clinDat)[grep(pattern = "Critical_Phenotype", x = row.names(clinDat), ignore.case = T)]

bloodDrawDates <- as.numeric(clinDat["N:M:CLIN:Data:Date_of_Blood_Collection_20150223__relative_to_Date_of_Birth",])
foldchange <- 0.5
pvalue <- 0.05

#####################################################################################################################

dayString <- c("AllDays","04Days","1Days")
today <- format(Sys.Date(), format="%Y%m%d")
dirs <- vector()
root <- paste0("/Users/chumphri/Projects/RNAseq/results/diffExpAnalysis/", today)
dir.create(path = root, showWarnings = T)

for(i in 1:length(dayString)){
  result <-  paste0(root, "/DE_blood_draw_", dayString[i])
  dir.create(path = result, showWarnings = T)
  dirs[i] <- result
}

for (day in 1:3) {
  covariates <- c("N:M:SURV:Data:Date_of_Birth__relative_to_Date_of_Birth", 
                  "C:M:ADMX:Data:Admix_80_Percent",
                  "N:M:CLIN:Data:Date_of_Blood_Collection_20150223__relative_to_Date_of_Birth")
  if (day == 1) {
    # then take them all 
    idx <- bloodDrawDates > -100000
  } 
  else if (day == 2) {
    idx <- bloodDrawDates >= 0 & bloodDrawDates <= 4
  
  } else if (day == 3) {
    idx <- bloodDrawDates == 1
  }
  idx[is.na(idx)] <- F
  clinMatFilt <- clinDat[,idx]  
  outdir <- dirs[day]; 
  
  for (target in targetPhenotypes) {
    probability<-find_probabilities(target = target, clin_Data = clinMatFilt)
    sample_size<-length(prob)
    for (repeatNumber in 1:10){
      a<-sample(names(e), replace=F, prob = probability)
      temp.mol<-molDat[ ,colnames(molDat) %in% names(a) ]
      targetString <- str_split(target, ":")[[1]][5]; #print(targetString)
      deTable <- diffExprFun(clinMat=clinMatFilt, dataMat=molDat,  targetPheno=target,
                           covarVec=covariates, FCThresh=foldchange, pValueThresh=pvalue,
                           writingDir=outdir, reps=10, cpus=1, robustFlag=F, writeTable=F)
    
     deTable$GeneID <- do.call("rbind", lapply(str_split(rownames(deTable), ":"), "[[", 5))
         
     numDE <- sum(deTable$adj.P.Val < pvalue & abs(deTable$logFC) >= foldchange)
     
     DEnumbers[repeatNumber] <- numDE
    
    print("*****************")
    print(paste0("Phenotype: ", targetString)) 
    print(paste0("Blood Draw Dates: ", dayString[day]))
    pos <- length(which(clinMatFilt[target,] == "true"))
    neg <- length(which(clinMatFilt[target,] == "false"))
    tot <-  ncol(clinMatFilt[target,])
    print(paste0("Number of samples [true]: ", pos))
    print(paste0("Number of samples [false]: ", neg))
    print(paste0("Total Number of samples: ", tot))
    prop <- pos/tot * 100
    print(paste0("Proportion of true samples: ", round(prop, digits = 1), "%" ))
    print(paste0("Number signficantly DE: ", sum(deTable$adj.P.Val < pvalue & abs(deTable$logFC) >= foldchange))) 
    print("*****************")
    }
    #foutstring <- add_date_tag(str_join(outdir, "/DE_RNA_", dayString[day], "_", targetString, "_val_targets"), ".txt")
    #write.table(deTable, quote=F, file=foutstring, sep="\t")
  }
  boxplot(x = DEnumbers, plot = T)
}