# miRNA Target Identification module using multiMiR Package Ru et al. 2014 
# Author: Nyasha Chambwe
# Date:   02/02/2015

require(stringr) || stop("Could not load package 'stringr'")

add_date_tag <- function(stringToTag, fileExtension){
  today <- Sys.Date()
  todayf <- format(today, format="%Y%m%d")
  return(paste(stringToTag, "_", todayf, sep = "", fileExtension))
}

parse_diff_exp_mirs <- function(setOfMiRs){
  mirIDs <- gsub(setOfMiRs, pattern = "N:M:MIRNA:mirID:", replacement = "")
  return(mirIDs)
}

# Load prequeried miRNA targets
validated_targets <- read.table("data/multiMir_validated_targets_20150213.txt")
# predicted_targets <- read.table("data/multiMir_predicted_targets_20150213.txt")

# for a given miRNA get validated targets from: mirecords, mirtarbase, and tarbase
get_validated_mir_targets <- function(mirOfInterest, targetPhenotypeName, writingDir){
  mirs <-parse_diff_exp_mirs(mirOfInterest)
  
  get_targets <- function (mirOfInterest) {
    #print(paste0("Searching for predicted gene targets of ", mirOfInterest))
    targetsOfInterest <- validated_targets[which(validated_targets$mature_mirna_id==mirOfInterest),]
    uniqueTargets <- unique(targetsOfInterest$target_symbol) 
    uniqueTargets <- droplevels(uniqueTargets)
    # print(paste0("Found ", length(uniqueTargets), " unique and validated targets for ", mirOfInterest))
    return(uniqueTargets)
  }
  
  listOfTargetGenes <- vector()
  for(mir in mirs){
    targets <- as.vector(get_targets(mir))
    listOfTargetGenes <- c(listOfTargetGenes, targets)
  }
  root <- add_date_tag(paste0("diffExp_miRNA_validated_target_genelist_", targetPhenotypeName), ".txt")
  filename <- paste0(writingDir, root)
  write.table(x = unique(listOfTargetGenes), file = filename, quote = F, row.names = F, col.names = F)
  return(listOfTargetGenes)  
}

# for a given miRNA get predicted targets from: TargetScan, Diana-microT, Elmmo, 
# Microcosm, Miranda, MirDB, Pictar, Pita
get_predicted_mir_targets <- function(mirOfInterest, targetPhenotypeName, writingDir){
  mirs <-parse_diff_exp_mirs(mirOfInterest)
  
  get_targets <- function (mirOfInterest) {
    #print(paste0("Searching for mRNA targets of ", mirOfInterest))
    targetsOfInterest <- predicted_targets[which(predicted_targets$mature_mirna_id==mirOfInterest),]
    uniqueTargets <- unique(targetsOfInterest$target_symbol) 
    uniqueTargets <- droplevels(uniqueTargets)
    #print(paste0("Found ", length(uniqueTargets), " unique and predicted targets for ", mirOfInterest))
    return(uniqueTargets)
  }
  
  listOfTargetGenes <- vector()
  for(mir in mirs){
    targets <- as.vector(get_targets(mir))
    listOfTargetGenes <- c(listOfTargetGenes, targets)
  }
  root <- add_date_tag(paste0("diffExp_miRNA_predicted_target_genelist_", targetPhenotypeName), ".txt")
  filename <- paste0(writingDir, root)
  write.table(x = unique(listOfTargetGenes), file = filename, quote = F, row.names = F, col.names = F)
  return(listOfTargetGenes)  
}