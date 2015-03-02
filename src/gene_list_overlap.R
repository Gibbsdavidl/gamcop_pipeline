require(GeneOverlap) || stop("Could not load package 'GeneOverlap'")
# list1 and list2 are the gene lists of interest: must be a character vector or a factor
# # additonal arguments are: 
#  genome.size - size of the gene universe - if NULL, then the total number of genes in that genome are used
#  spec - species  can be one of mm9.gene, hg19.gene, rn4.gene
get_list_overlap_significance <- function(list1, list2, geneUniverse, species){
  go.obj <- newGeneOverlap(listA=list1, listB = list2, genome.size=geneUniverse, spec=species)
  go.obj <- testGeneOverlap(go.obj)
  return(getPval(go.obj))
}

get_list_overlap_size <- function(list1, list2, geneUniverse, species){
  go.obj <- newGeneOverlap(listA=list1, listB = list2, genome.size=geneUniverse, spec=species)
  go.obj <- testGeneOverlap(go.obj)
  return(length(getIntersection(go.obj)))
}

hypergeomTest <- function(genelist, otherlist, universe) {
  # Each variable needs to be a character vector.
  stopifnot(is.character(genelist), is.character(otherlist), is.character(universe))
  #phyper(x,m,n,k)
  #x, vector of quantiles representing the number of white balls drawn without replacement from an urn which contains both black and white balls.
  #m   the number of white balls in the urn.
  #n   the number of black balls in the urn.
  #k	 the number of balls drawn from the urn.
  
  x <- length(intersect(genelist, otherlist))  # sig and in pathway for example
  
  # Source: http://stats.stackexchange.com/questions/16247/calculating-the-probability-of-gene-list-overlap-between-an-rna-seq-and-a-chip-c
  # This would be the enrichment p val
  pval <- phyper(x-1, length(otherlist), (length(universe) - length(otherlist)), length(genelist), lower.tail = F)
  
  # Numerical parameters in order:
  # (success-in-sample, success-in-bkgd, failure-in-bkgd, sample-size).
  return(pval)
}

fisherTestEnrichment <- function(list1, list2, geneUniverse){
  a <- length(intersect(list1, list2))
  b <- length(setdiff(list1, list2))
  c <- length(setdiff(list2, list1))
  d <- length(geneUniverse) - length(union(list1, list2))
  cont.table <- matrix(c(a,b,c,d), nrow=2)
  fisherP.htest <- fisher.test(x = cont.table, alternative = "greater")
  return(fisherP.htest$p.value)
}