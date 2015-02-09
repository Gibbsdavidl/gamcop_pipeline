# miRNA Target Identification module using multiMiR Package
# Ru et al. 2014 
# Author: Nyasha Chambwe
# Date:   02/02/2015
require(XML) || stop("Could not load package 'XML'")
require(RCurl) || stop("Could not load package 'RCurl'")
require(multiMiR) || stop("Could not load package 'multiMiR'")
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

# for a given miRNA get validated targets
# Searches mirecords, mirtarbase, and tarbase for hits
get_validated_mir_targets <- function(mirOfInterest, targetPhenotypeName, writingDir){
  mirs <-parse_diff_exp_mirs(mirOfInterest)
  
  get_targets <- function (mirOfInterest) {
    print(paste0("Searching for mRNA targets of ", mirOfInterest))
    table <- get.multimir(mirna=mirOfInterest, summary=TRUE)
    uniqueTargets <- unique(table$validated$target_symbol) 
    print(paste0("Found ", length(uniqueTargets), " unique and validated targets for", mirOfInterest))
    return(uniqueTargets)
  }
  
  listOfTargetGenes <- sapply(mirs, function(x) get_targets(x))
  setOfTargetGenes <- unlist(x =listOfTargetGenes, use.names = F)
  
  root <- add_date_tag(paste0("diffExp_miRNA_target_genelist_", targetPhenotypeName), ".txt")
  filename <- paste0(writingDir, root)
  write.table(x = setOfTargetGenes, file = filename, quote = F, row.names = F, col.names = F)
  
}