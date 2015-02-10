# Module to compare the overlap of differentially expressed/methylated genes either with each other or with curated gene sets
# Author: Nyasha Chambwe
# Date: 02/06/2015
source("../../../../src/draw_venn_diagram.R")
source("../../../../src/gene_list_overlap.R")

add_date_tag <- function(stringToTag, fileExtension){
  today <- Sys.Date()
  todayf <- format(today, format="%Y%m%d")
  return(paste(stringToTag, "_", todayf, sep = "", fileExtension))
}

# Arguements
# targetPhenotype: Description of Target Phenotype for which genes are differentially expressed
# rna: List of Differentially Expressed Genes : characterVector
# mir: List of Differnetially Expressed miRs : characterVector
# meth: List of Differnetially Methylated Genes
# geneList: Any other list (dbPTB gene list)
# writingDir: directory to write the results 
get_de_genes_overlap <- function(rna, mir, meth, geneList, targetPhenotypeName, writingDir){
  
  perform_overlap_analysis <- function(listOfItemsToCompare, overlapStats, targetPhenotypeName) {
    comparisonName <- paste(names(listOfItemsToCompare)[1], names(listOfItemsToCompare)[2], sep="_")
    inCommon <- get_list_overlap_size(listOfItemsToCompare[[1]], listOfItemsToCompare[[2]], NULL, species="hg19.gene")
    pval <- sprintf("%.2e", get_list_overlap_significance(listOfItemsToCompare[[1]], listOfItemsToCompare[[2]], 
                                                          NULL, species="hg19.gene"))
    overlapStats[nrow(overlapStats)+1,] <- c(comparisonName, length(listOfItemsToCompare[[1]]), length(listOfItemsToCompare[[2]]), inCommon, pval)
    root <- paste("overlap_venn", targetPhenotypeName, names(listOfItemsToCompare)[1], names(listOfItemsToCompare)[2], sep = "_")
    filename <- add_date_tag(root, ".tiff")
    filename_path <- paste0(writingDir, filename)
    draw_venn_diagram(listOfItemsToCompare, filename_path)
    return(overlapStats)
  }
  
  resultsTable <- data.frame( "Comparison"=character(), "List 1 Size" = integer(), "List 2 Size" = integer(), 
                              "Overlap Size" = integer(), "P-value"=integer(), stringsAsFactors=FALSE)
  
  # pairwise comparisons between molecular data types
  if(!is.null(rna) & !is.null(mir)){
  compList <- list("mRNA_Expression"=rna, "miRNA_Targets"=mir )
  resultsTable <- perform_overlap_analysis(compList, resultsTable, targetPhenotypeName = targetPhenotypeName)
  }
  
  if(!is.null(rna) & !is.null(meth)){
  compList <- list("mRNA_Expression"=rna, "Methylation"=meth )
  resultsTable <- perform_overlap_analysis(compList, resultsTable, targetPhenotypeName)
  }
  
  if(!is.null(meth) & !is.null(mir)){
  compList <- list("Methylation"=meth, "miRNA_Targets"=mir )
  resultsTable <- perform_overlap_analysis(compList, resultsTable, targetPhenotypeName)
  }
  
  # pairwise comparisons between molecular data types and PTB Gene List
  compList <- list("mRNA_Expression"=rna, "GeneList"= geneList )
  resultsTable <- perform_overlap_analysis(compList, resultsTable, targetPhenotypeName)
  
  compList <- list("miRNA_Targets"=mir, "GeneList"= geneList)
  resultsTable <- perform_overlap_analysis(compList, resultsTable, targetPhenotypeName)
  
  if(!is.null(meth)){
  compList <- list("Methylation"=meth, "GeneList"= geneList )
  resultsTable <- perform_overlap_analysis(compList, resultsTable, targetPhenotypeName)
  }
  
  return(resultsTable)
}
