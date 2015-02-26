# Make sure gene symbols are approved HGNC symbols
# Author: Nyasha Chambwe
# Date:   20150226
require(HGNChelper) || stop("Could not load package 'GeneOverlap'")

# Input: character vector of gene names 
# Returns character vector (same length as input vector) with HGNC approved gene names
# If no approved gene name is found - <NA> is returned
getHgncSymbol <- function(geneSymbolsToCheck){
  dat <- checkGeneSymbols(geneSymbolsToCheck, unmapped.as.na=TRUE)
  dat$Suggested.Symbol
}

# Returns how many genes without matched Aliases
getCountNonHGNC <- function(geneSymbolsToCheck){
  dat <- checkGeneSymbols(geneSymbolsToCheck, unmapped.as.na=TRUE)
  print(dat[is.na(dat$Suggested.Symbol),])
  return(sum(is.na(dat$Suggested.Symbol)))
}