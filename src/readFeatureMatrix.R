


# For example
# fmfile <- "2015_01_14_genomic_clinical_hilevel.fm"
# dataMat <- readFeatureMatrix(methfile)
# and you get a matrix with *real* column and row names.


readFeatureMatrix <- function(pathToFM=NA)
{
  require("stringr")
  MatA <- read.delim(pathToFM, header=FALSE, stringsAsFactors=FALSE)
  MatB <- MatA[-1,-1]
  rownames(MatB) <- MatA[-1,1]
  colnames(MatB) <- MatA[1,-1]
  as.data.frame(MatB)
  dataType <- str_sub(rownames(MatB), 1, 1)
  if (all(dataType == "N")) {
      rownamesMatB <- rownames(MatB)
      MatB <- apply(MatB, MARGIN=2, FUN=function(a) as.numeric(a))
      rownames(MatB) <- rownamesMatB
  }
  MatB
}
