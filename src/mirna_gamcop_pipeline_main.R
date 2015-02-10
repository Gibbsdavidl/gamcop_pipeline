# Script to run miRNA-Seq differential expression analysis using the unified GAMCOP Pipeline
# Date: 02/02/2015
# Author: Nyasha Chambwe

# Source Functions and R Packages -----------------------------------------
# Source functions from: https://github.com/Gibbsdavidl/gamcop_pipeline.git
# Developed in collaboration between Crystal Humphries, Dave Gibbs and Nyasha Chambwe
source("src/diffExprFun.R")
source("src/differential.R")
source("src/readFeatureMatrix.R")
source("src/functions_visualize_diff_exp.R")
source("src/visualize_diff_exp.R")
source("src/additional_data_parsing.R")
source("src/identify_mir_targets.R")

# Returns limma diffExp table of signficant hits according to pval and fc thresholds
# Using names may contribute to misaligned headers
filter_diffExp <- function(pathToDeTable, FCThresh, pValueThresh){
  print(paste0("Reading in file: ", pathToDeTable))
  tab <- read.table(pathToDeTable, header = T, sep="")
  stopifnot( "logFC" %in% names(tab), "adj.P.Val" %in% names(tab))
  diffExpTab <- tab[abs(tab$logFC) >=FCThresh & tab$adj.P.Val <=pValueThresh,]
  return(diffExpTab)
}
# Load Input Files --------------------------------------------------------

# Load Input Files
# Clinical data file
clinDat <-  read.table("/Volumes/~nchambwe/inova/projects/gamcop/smrna-analysis/data/clinical/data_CLIN_20150203.fm", 
                       check.names=F, header=T, row.names = 1)
molDat <- read.table("/Volumes/~nchambwe/inova/projects/gamcop/smrna-analysis/data/molecular/small_rna_seq/Study_101/batch_summaries/data_MIRSeq_LIMMA_corrected_20141201.fm", 
                     header=T, sep="\t", row.names=1, check.names = F)

targetPhenotype <- "N:NB:CLIN:Critical_Phenotype:TermCategory"

#target phenotypes (all twelve are here)
targetPhenotypes <- c("B:NB:CLIN:Critical_Phenotype:1n2v4",
                      "B:NB:CLIN:Critical_Phenotype:1v4",
                      "N:NB:CLIN:Critical_Phenotype:Gestational_Age_at_Delivery",
                      "B:NB:MRGE:Critical_Phenotype:History_of_Preterm_Birth",
                      "B:NB:MRGE:Critical_Phenotype:Hypertension_Related",
                      "B:NB:MRGE:Critical_Phenotype:Immune_Related",
                      "B:M:CLIN:Critical_Phenotype:Incompetent_Shortened_Cervix",
                      "B:M:CLIN:Critical_Phenotype:Inova_Idiopathic_NA",
                      "B:NB:MRGE:Critical_Phenotype:Placenta_Related",
                      # "B:M:CLIN:Critical_Phenotype:Preeclampsia/Eclampsia",
                      "B:NB:CLIN:Critical_Phenotype:Preterm",
                      "B:NB:MRGE:Critical_Phenotype:Prom_Related",
                      "N:NB:CLIN:Critical_Phenotype:TermCategory",
                      "B:NB:MRGE:Critical_Phenotype:Uterine_Related"
)

covariates = c("N:F:SURV:Data:Date_of_Birth__relative_to_Date_of_Birth", "C:F:ADMX:Data:Admix_80_Percent") 

# Differential Expression Analysis: All Samples ---------------------------
foldChange <- 0.25
pvalue <- 0.05
resultsLocation <- "/Volumes/~nchambwe/inova/projects/gamcop/gamcop_pipeline/results/diffExpAnalysis/20150209/all_samples/"
gamcop<- targetPhenotypes

sapply(gamcop, function(x) diffExprFun(clinMat = clinDat, dataMat=molDat, targetPheno=x, 
                                       covarVec=covariates, FCThresh=foldChange, pValueThresh=pvalue, 
                                       writingDir=resultsLocation))

diffExpTables <- list.files(path = resultsLocation, pattern="DE_MIRNA_*", all.files = F, recursive = T)

for(i in 1:length(gamcop)){
  sigMirsTable <- filter_diffExp(pathToDeTable =paste0(resultsLocation, diffExpTables[i]), FCThresh = foldChange, pValueThresh = pvalue)
  mirs <- row.names(sigMirsTable)
  comparison <- unlist(str_split(gamcop[i], pattern = ":"))[[5]]
  print(paste0("The number of differentially expressed mirs for ", comparison, " is: ", length(mirs)))
  
  if(length(mirs)>=1){
    get_validated_mir_targets(mirOfInterest = mirs, targetPhenotypeName =comparison , 
                              writingDir =resultsLocation)
  }
  print("Completed: miRNA Target identification")
}

# Differential Expression Analysis: Around Labor Day (Day 0) to (Day 4) of Sample Collection --------------

#clin.rowname = rownames(Clin.dat)['desired trait']
#traits        = vector of matching items wanted (e.g. blood draw date ranges from -31  to 309 days. To get days 1 and 2, 
bloodDrawDates <- "N:M:CLIN:Data:Date_of_Blood_Collection__relative_to_Date_of_Birth"
molDat_blood_draw_0_4 <- data_subset_Clin(Mol.dat.frame = molDat, Clin.dat.frame = clinDat, clin.rowname =bloodDrawDates , traits = c(0:4))

resultsLocation <- "/Volumes/~nchambwe/inova/projects/gamcop/gamcop_pipeline/results/diffExpAnalysis/20150209/blood_draw_0_4/"

sapply(gamcop, function(x) diffExprFun(clinMat = clinDat, dataMat=molDat_blood_draw_0_4, targetPheno=x, 
                                       covarVec=covariates, FCThresh=foldChange, pValueThresh=pvalue, 
                                       writingDir=resultsLocation))

diffExpTables <- list.files(path = resultsLocation, pattern="DE_MIRNA_*", all.files = F, recursive = T)

for(i in 1:length(gamcop)){
  sigMirsTable <- filter_diffExp(pathToDeTable =paste0(resultsLocation, diffExpTables[i]), FCThresh = foldChange, pValueThresh = pvalue)
  mirs <- row.names(sigMirsTable)
  comparison <- unlist(str_split(gamcop[i], pattern = ":"))[[5]]
  print(paste0("The number of differentially expressed mirs for ", comparison, " is: ", length(mirs)))
  
  if(length(mirs)>=1){
    get_validated_mir_targets(mirOfInterest = mirs, targetPhenotypeName =comparison , 
                              writingDir =resultsLocation)
  }
  print("Completed: miRNA Target identification")
}
