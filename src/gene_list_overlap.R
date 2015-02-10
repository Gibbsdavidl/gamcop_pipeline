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
