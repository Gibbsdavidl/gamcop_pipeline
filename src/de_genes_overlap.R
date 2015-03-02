# Module to compare the overlap of differentially expressed/methylated genes either with each other or with curated gene sets
# Author: Nyasha Chambwe
# Date: 02/06/2015
source("../../../../src/draw_venn_diagram.R")
source("../../../../src/gene_list_overlap.R")

formatPval <- function(num){
  sprintf("%.2e", num)
}

add_date_tag <- function(stringToTag, fileExtension){
  today <- Sys.Date()
  todayf <- format(today, format="%Y%m%d")
  return(paste(stringToTag, "_", todayf, sep = "", fileExtension))
}

write_overlap_venn <- function (targetPhenotypeName, listOfItemsToCompare, writingDir) {
  root <- paste("overlap_venn", targetPhenotypeName, names(listOfItemsToCompare)[1], names(listOfItemsToCompare)[2], sep = "_")
  filename <- add_date_tag(root, ".tiff")
  filename_path <- paste0(writingDir, filename)
  draw_venn_diagram(listOfItemsToCompare, filename_path)
}

perform_overlap_signficance_analysis <- function(listOfItemsToCompare, overlapStats, targetPhenotypeName, universe) {
  list1.size <- length(listOfItemsToCompare[[1]])
  list2.size <- length(listOfItemsToCompare[[2]])
  universe.size <- length(universe)
  expectedOverlap <- round((list1.size*list2.size)/universe.size)
  inCommon <- get_list_overlap_size(listOfItemsToCompare[[1]], listOfItemsToCompare[[2]], universe.size, species="hg19.gene")
  goPval <- get_list_overlap_significance(listOfItemsToCompare[[1]], listOfItemsToCompare[[2]], universe.size, species="hg19.gene")
  hyperP <- hypergeomTest(as.character(listOfItemsToCompare[[1]]), as.character(listOfItemsToCompare[[2]]), as.character(universe))
  fisherP_manual <-fisherTestEnrichment(listOfItemsToCompare[[1]], listOfItemsToCompare[[2]], as.character(universe))
  overlapStats[nrow(overlapStats)+1,] <- c(names(listOfItemsToCompare)[1], list1.size, 
                                           names(listOfItemsToCompare)[2], list2.size, 
                                           universe.size, expectedOverlap, inCommon, 
                                           formatPval(goPval), formatPval(goPval), formatPval(hyperP))
  return(overlapStats)
}


# Arguements
# targetPhenotype: Description of Target Phenotype for which genes are differentially expressed
# rna: List of Differentially Expressed Genes : characterVector
# mir: List of Differnetially Expressed miRs : characterVector
# meth: List of Differnetially Methylated Genes
# geneList: Any other list (dbPTB gene list)
# writingDir: directory to write the results 

get_de_genes_overlap <- function(rna, mir, meth, geneList, targetPhenotypeName, writingDir, rnaUniverse, 
                                 mirUniverse, methUniverse, writeVenn=FALSE){
  
  
  resultsTable <- data.frame( "GeneList1"=character(), "List1 " = integer(), "GeneList2"=character(), "List2" = integer(), 
                              "Universe"= integer(), "Expected" = integer(), "Actual" = integer(), "Pval_F"=numeric(), "Pval_O"=numeric(), "HyperP"=numeric(), stringsAsFactors=FALSE)
  
  # pairwise comparisons between molecular data types
  if(!is.null(rna) & !is.null(mir)){
    universeGeneList <- intersect(rnaUniverse, mirUniverse)
    rnaUni <- rna[rna %in% universeGeneList]
    mirUni <- mir[mir %in% universeGeneList]
    compList <- list("mRNA_Expression"=rnaUni, "miRNA_Targets"=mirUni )
    resultsTable <- perform_overlap_signficance_analysis(compList, resultsTable, targetPhenotypeName = targetPhenotypeName, universeGeneList)
    if(writeVenn){write_overlap_venn(targetPhenotypeName, compList, writingDir)}
  }    
  
  if(!is.null(rna) & !is.null(meth)){
    universeGeneList <- intersect(rnaUniverse, methUniverse)
    rnaUni <- rna[rna %in% universeGeneList]
    methUni <- meth[meth %in% universeGeneList]
    compList <- list("mRNA_Expression"=rnaUni, "Methylation"=methUni)
    resultsTable <- perform_overlap_signficance_analysis(compList, resultsTable, targetPhenotypeName, universeGeneList)
    if(writeVenn){write_overlap_venn(targetPhenotypeName, compList, writingDir)}
  }
  
  if(!is.null(meth) & !is.null(mir)){
    universeGeneList <- intersect(methUniverse, mirUniverse)
    universe <- length(universeGeneList)
    methUni <- meth[meth %in% universeGeneList]
    mirUni <- mir[mir %in% universeGeneList]
    compList <- list("Methylation"=methUni, "miRNA_Targets"=mirUni)
    resultsTable <- perform_overlap_signficance_analysis(compList, resultsTable, targetPhenotypeName, universeGeneList)
    if(writeVenn){write_overlap_venn(targetPhenotypeName, compList, writingDir)}
  }
  
  # pairwise comparisons between molecular data types and PTB Gene List
  
  if(!is.null(rna) & !is.null(geneList)){
    #universeGeneList <- intersect(rnaUniverse, geneList)
    universeGeneList <- rnaUniverse
    rnaUni <- rna[rna %in% universeGeneList]
    geneListUni <- geneList[geneList %in% universeGeneList]
    compList <- list("mRNA_Expression"=rnaUni, "dbPTBGeneList"= geneListUni )
    resultsTable <- perform_overlap_signficance_analysis(compList, resultsTable, targetPhenotypeName, universeGeneList)
    if(writeVenn){write_overlap_venn(targetPhenotypeName, compList, writingDir)}
  }
  
  if(!is.null(mir) & !is.null(geneList)){
    #universeGeneList <- intersect(mirUniverse, geneList)
    universeGeneList <- mirUniverse
    universe <- length(universeGeneList)
    mirUni <- mir[mir %in% universeGeneList]
    geneListUni <- geneList[geneList %in% universeGeneList]
    compList <- list("miRNA_Targets"=mirUni, "dbPTBGeneList"= geneListUni)
    resultsTable <- perform_overlap_signficance_analysis(compList, resultsTable, targetPhenotypeName, universeGeneList)
    if(writeVenn){write_overlap_venn(targetPhenotypeName, compList, writingDir)}
  }
  
  if(!is.null(meth) & !is.null(geneList)){
    #universeGeneList <- intersect(methUniverse, geneList)
    universeGeneList <- methUniverse
    universe <- length(universeGeneList)
    methUni <- meth[meth %in% universeGeneList]
    geneListUni <- geneList[geneList %in% universeGeneList]
    compList <- list("Methylation"=methUni, "dbPTBGeneList"= geneListUni )
    resultsTable <- perform_overlap_signficance_analysis(compList, resultsTable, targetPhenotypeName, universeGeneList)
    if(writeVenn){write_overlap_venn(targetPhenotypeName, compList, writingDir)}
  }
  
  filename <- add_date_tag(stringToTag = paste0("overlap_table_", targetPhenotypeName), fileExtension = ".txt")
  write.table(resultsTable, file = filename, quote = F, col.names = T, row.names = F)
  return(resultsTable)
}