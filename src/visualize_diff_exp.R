# Module to visualize the results of differential expression analysis
# Author: Nyasha Chambwe
# Date: 20150130
source("src/functions_visualize_diff_exp.R")
require(ggplot2) || stop("Could not load package 'ggplot2'")
require(reshape2)|| stop("Could not load package 'reshape2'")
require(stringr) || stop("Could not load package 'stringr'")
# very annoying everytime I try to merge df1, df2 by row names have to constantly replace the rownames
# merge df1 and df2 and keep the rownames the same
# only produce output for those rows that could be merged
nyMerge <- function(df1, df2){
  #stopifnot(is.data.frame(df1), is.data.frame(df2))
  merged_df <- merge(df1, df2, by=0)
  row.names(merged_df) <- merged_df$Row.names
  merged_df <- merged_df[,-1]
  return(merged_df)
}



# Inputs
# topTable - the full table that results from differential expression analysis
# topK - how many top genes to make individual panels for
# targetPheno - phenotype for comparison
# FCThresh  - fold change threshold 
# inputs passed from main function: dataMat, targetPheno, FCThresh, pValueThresh, writingDir
visualize_diff_exp <- function(clinMat, dataMat, topTable, topk=1, targetPheno, FCThresh, pValueThresh, writingDir){
  # Human readable phenotype name
  phenotypeName <- unlist(str_split(targetPheno, ":"))[[5]]
  diffExpTable <- topTable[topTable$adj.P.Val<=pValueThresh & abs(topTable$logFC>=FCThresh),]
  diffExpGenes <- row.names(diffExpTable)
  numDiffExp <- length(diffExpGenes)
  print(paste("The number of differentially expressed genes is: ", numDiffExp))
  # Only plot things if you actually have significant differences
  if(numDiffExp > 0){
    dataSource <- unlist(str_split(diffExpGenes[[1]], ":"))[[3]]
    # N/B/C?
    phenotypeDataType <- unlist(str_split(targetPheno, ":"))[[1]]
    
    # Volcano plots
    make_volcano_plot(topTable, dataSource, pValueThresh, FCThresh, phenotypeName)
    
    # Make boxplots of the top k expressed genes (defaults to 5)
    minTopFeatures <- min(numDiffExp, topk)
    
    # Create ggplot2 input data
    for(gene in diffExpGenes[1:minTopFeatures]){
      geneOfInterest <- gene
      geneExpressionLevel <- as.data.frame(t(dataMat[geneOfInterest,]))
      classTable <- as.data.frame(clinMat[,targetPheno], row.names = row.names(clinMat))
      dat.obj <- nyMerge(classTable, geneExpressionLevel)
      expression_by_phenotype_boxplots(dat.obj, unlist(str_split(geneOfInterest, ":"))[[5]], phenotypeName, dataSource, writingDir)
    }
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