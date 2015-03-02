
dataTyper <- function(type, dat) 
{
  if (type == "N") {
    return(as.numeric(dat))
  } else if (type == "B") {
    return(as.numeric(as.factor(dat)))
  } else if (type == "C") {
    return(as.factor(dat))
  }
  return(dat)
}

diffExprErrorCheck <- function(clinMat, dataMat, targetPheno=NA, covarVec=NA, 
                               FCThresh=NA, pValueThresh=NA, writingDir=".")
{
  require(stringr)
  if (is.na(targetPheno)) {
    print("diffExprFun Error: Please specify the target phenotype.")
    return(NA)
  }
  if (!is.numeric(FCThresh)) {
    print("diffExprFun Error: FCThresh must be numeric")
    return(NA)
  }
  if ((!is.numeric(pValueThresh)) | pValueThresh < 0 | pValueThresh > 1) {
    print("diffExprFun Error: pValue must be between 0 and 1")
    return(NA)
  }
  if (is.na(targetPheno)) {
    print("diffExprFun Error: Please specify the target phenotype.")
    return(NA)
  }
  if(file.info(writingDir)[1,"isdir"] != TRUE) {
    print("diffExprFun Error: This is a directory. Please specify a file name.")
    return(NA)
  }
  ms <- colnames(dataMat)
  ms <- str_sub(ms, 1,7)     ########## compare the first 7 chars "101-474" ##########
  idx <- c(covarVec, targetPheno)
  ns <- str_sub(colnames(clinMat), 1,7)

  if(sum(idx %in% rownames(clinMat)) < length(idx)) {
      print("diffExprFun Error: Row names not found in clinical matrix.")
      return(NA)
  }
  
  clinMatDataFrame <- t(clinMat[idx, ns %in% ms])
  
  #if (nrow(clinMatDataFrame) != length(ms)) {
  #  print("diffExprFun Error: nrow of clinMatSubset doesn't equal the number of samples.")
  #  return()
  #}
    
  dataTypes <- str_sub(idx, 1,1)
  clinMatTyped <- as.data.frame(sapply(1:length(dataTypes), FUN=function(a) dataTyper(dataTypes[a], clinMatDataFrame[,a])))
  colnames(clinMatTyped) <- colnames(clinMatDataFrame)
  row.names(clinMatTyped) <- row.names(clinMatDataFrame)
  if(all(dim(clinMatTyped) != c(length(ms), length(dataTypes)))) {
    print("diffExprFun Error: clinMat formatting failed somehow..")
    return(NA)
  }
  
  design <- model.matrix(~., data=clinMatTyped)    
  # design
  
  itemsToReturn <- list(designMat=design, clinDF=clinMatTyped)
  itemsToReturn
}

diffExprFun <- function(clinMat, dataMat, targetPheno=NA, covarVec=NA, FCThresh=NA, pValueThresh=NA, writingDir="./") 
{
  # Assumes that:
  # clinMat and dataMat have the same column names, i.e. "101-479-M" .. R data.frame column names in are in place.
  # clinMat and datMat have rownames .... preprocessing will make sure the rowname column is eliminated.
  # the target phenotype is a character string like "B:CLIN:Preterm:NB::::" 
  # the covarVec is a vector of strings matching 1st-column-names in clinMat
  # FCthresh is a number like 1.2
  # pValueThresh is a number between 0 and 1
  # the tables will be written to the current directory, and returned to pass on in the pipeline
  
  # Check for Input Errors --------------------------------------------------
  errorCheckOutputs <- diffExprErrorCheck(clinMat, dataMat, targetPheno, covarVec, FCThresh, pValueThresh, writingDir)
  # design matrix
  designTable <- errorCheckOutputs[[1]]
  # clinical data matrix (samples with molecular data)
  clinMatFiltered <- errorCheckOutputs[[2]]
  if (is.na(designTable)) {
    # error!
    return(NA)
  }
  print("Completed: Input error check ")
  # Differential expression analysis ----------------------------------------
  
  topTable <- differential(designTable, dataMat, covarVec, writingDir)
  print("Completed: Differential expression testing ")
  
  # Visualize results of differential expression analysis -------------------
  
  visualize_diff_exp(clinMatFiltered, dataMat, topTable, topk=5, targetPheno, FCThresh, pValueThresh, writingDir)
    
  print("Completed: Visualization of differential expression test results ")
  return(topTable)
}


bootDiffFun <- function(clinMat, dataMat, targetPheno=NA, covarVec=NA, 
                        FCThresh=NA, pValueThresh=NA, writingDir="./", 
                        reps=100, cpus=2, writeTable=F)
{
  # Assumes that:
  # clinMat and dataMat have the same column names, i.e. "101-479-M" .. R data.frame column names in are in place.
  # clinMat and datMat have rownames .... preprocessing will make sure the rowname column is eliminated.
  # the target phenotype is a character string like "B:CLIN:Preterm:NB::::" 
  # the covarVec is a vector of strings matching 1st-column-names in clinMat
  # FCthresh is a number like 1.2
  # pValueThresh is a number between 0 and 1
  # the tables will be written to the current directory, and returned to pass on in the pipeline
  
  # Check for Input Errors --------------------------------------------------
  errorCheckOutputs <- diffExprErrorCheck(clinMat, dataMat, targetPheno, covarVec, FCThresh, pValueThresh, writingDir)
  # design matrix
  designTable <- errorCheckOutputs[[1]]
  # clinical data matrix (samples with molecular data)
  clinMatFiltered <- errorCheckOutputs[[2]]
  if (is.na(designTable)) {
    # error!
    return(NA)
  }
  print("Completed: Input error check ")
  # Differential expression analysis ----------------------------------------
  
#  topTable <- differential(design, dataMat, covarVec, writingDir)
  topTable <- bootDiff(designTable, dataMat, covarVec, writingDir, 5, reps, cpus, writeTable)
  print("Completed: Differential expression testing ")
  
  # Visualize results of differential expression analysis -------------------
  
  visualize_diff_exp(clinMatFiltered, dataMat, topTable, topk=5, targetPheno, FCThresh, pValueThresh, writingDir)
  
  print("Completed: Visualization of differential expression test results ")
  return(topTable)
}






