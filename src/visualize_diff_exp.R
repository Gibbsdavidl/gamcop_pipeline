# Module to visualize the results of differential expression analysis
# Author: Nyasha Chambwe
# Date: 20150130
source("functions_visualize_diff_exp.R")
require(ggplot2) || stop("Could not load package 'ggplot2'")
require(reshape2)|| stop("Could not load package 'reshape2'")

# Inputs
# topTable - the full table that results from differential expression analysis
# topK - how many top genes to make individual panels for
# targetPheno - phenotype for comparison
# FCThresh  - fold change threshold 
# inputs passed from main function: dataMat, targetPheno, FCThresh, pValueThresh, writingDir
visualize_diff_exp <- function(itemsToReturn, dataMat, topTable, topk=5, targetPheno, FCThresh, pValueThresh, writingDir){
  clinMat <- itemsToReturn[[2]]
  phenotypeOfInterest <- unlist(str_split(targetPheno, ":"))[[5]]
  dataType <- unlist(str_split(targetPheno, ":"))[[1]]
  dataSource <- unlist(str_split(targetPheno, ":"))[[1]]
  diffExpTable <- topTable[topTable$adj.P.Val<=pValueThresh & abs(topTable$logFC>=FCThresh),]
  diffExpGenes <- row.names(diffExpTable)
  numDiffExp <- length(diffExpGenes)
  
  # Only plot things if you actually have differences
  if(numDiffExp > 0){
    # Volcano plots
    make_volcano_plot(topTable, dataSource, pValueThresh, FCThresh, phenotypeOfInterest)
    
    # Make boxplots of the top k expressed genes (defaults to 5)
    minTopFeatures <- min(numDiffExp, topk)
    
    # Create ggplot2 input data
    
    #geneOfInterest <- diffExpGenes[1]
    #geneExpressionLevel <- t(dataMat[geneOfInterest,])
    #class_table <- t(clinMat[phenotypeOfInterest,])
    #dat.obj <- nyMerge(class_table, expression_of_gene_of_interest)
    
    # @TODO 
    # Plot heatmap
    # colors change depending on data source (paper convention)
    # if(dataSource=="METH", dataSource=="RNASEQ", dataSource=="MIRNA"
    # print(paste(dataSource, "is not a valid data type for this analysis!", sep = " "))
  } else {
    print("There are no signficant differentially expressed/methylated genes for the selected P value 
          and fold change thresholds") 
  }
}
