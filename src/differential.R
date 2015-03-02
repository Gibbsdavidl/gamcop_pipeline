differential <- function(designTable, dataMat, covarVec, writingDir, FCThresh, pValueThresh,...){
   require(limma)
  
  FixedDataMatrix<-dataMat[, colnames(dataMat) %in% rownames(designTable)]
  
  #create file name for data output
  DataType       = strsplit(rownames(FixedDataMatrix)[1], split=":")[[1]][3]
  targetPheno <- unlist(str_split(colnames(designTable)[ncol(designTable)], ":"))[[5]]  # last col of designTable
  targetPheno <- str_replace_all(string = targetPheno, pattern = "`", replacement = "")
  targetPheno <- gsub("/","_", targetPheno)
  print(targetPheno)
  OutputFile     = paste0("/DE_", DataType, "_", targetPheno)
  OutputFileName = add_date_tag(OutputFile, ".txt")
  OutputFile_dir = paste0(writingDir, OutputFileName)
  
  # make sure that the samples are the same in both the dataMatrix and designTable matrix
  CheckData<-table(colnames(FixedDataMatrix) %in% rownames(designTable))
  
  if (length(CheckData) == 1 & names(CheckData)[1] == "TRUE"){
    
    # run limma
    fit1<-lmFit(FixedDataMatrix, designTable)
    fit1<-eBayes(fit1)
    
    #gather relevant covariates
    cov_length = length(covarVec)+2
    DesignVariables  = dim(designTable)[2]
    Coef_of_Interest = paste(cov_length,DesignVariables, sep=":")
    
    #Write and Return Table of DE/M variables and the corresponding statistics
    print(cov_length:DesignVariables)
    Complete_table<-topTable(fit1, n=Inf, coef= cov_length:DesignVariables)
    ID_name<-paste0("N:M:",DataType, ":Data:")
    Complete_table$GeneName <- gsub(ID_name, "", rownames(Complete_table))
    write.table(Complete_table, file=OutputFile_dir, quote = F, row.names = T, sep="\t")
    return(Complete_table)
  }
  
  else{
      stop("Sample IDs in Design Matrix do not match Sample IDs in Data Matrix")
  }
}


runLimma <- function(DataMatrix, designTable, covarVec)
{  
  fit1<-lmFit(DataMatrix, designTable)
  fit1<-eBayes(fit1)
  
  #gather relevant covariates
  cov_length = length(covarVec)+2
  DesignVariables  = dim(designTable)[2]
  Coef_of_Interest = paste(cov_length, DesignVariables, sep=":")
  
  #Write and Return Table of DE/M variables and the corresponding statistics
  print(cov_length:DesignVariables)
  Complete_table<-topTable(fit1, n=Inf, coef= cov_length:DesignVariables)
  Complete_table
}


bootLimmaStat <- function(FixedDataMatrix, CompleteTable, designTable, covarVec, j, repidx)
{
  # i == which rep it is
  # j == index of : {logFC    AveExpr        t      P.Value    adj.P.Val        B}    
  ids <- intersect(colnames(FixedDataMatrix), rownames(designTable))
  FixedDataMatrix <- FixedDataMatrix[,colnames(FixedDataMatrix) %in% ids]
  designTable <- designTable[rownames(designTable) %in% ids, ]
  N <- ncol(FixedDataMatrix)
  X <- mat.or.vec(nr = nrow(FixedDataMatrix), nc=length(repidx))
  for (i in 1:length(repidx)) {
    sampleIdx <- sample(1:N, size=N, replace=T)
    SampledDataMatrix    <- FixedDataMatrix[,sampleIdx]
    SampledDesign        <- designTable[sampleIdx,]
    SampledCovar         <- 
    SampledCompleteTable <- runLimma(SampledDataMatrix, SampledDesign, covarVec)
    SampledCompleteTable <- SampledCompleteTable[match(table = rownames(SampledCompleteTable), 
                                                     x = rownames(CompleteTable)),]
    X[,i] <- SampledCompleteTable[,j]
  }
  rownames(X) <- rownames(CompleteTable)
  X
}


runBootstrap <- function(FixedDataMatrix, CompleteTable, designTable, covarVec, cpus, j, reps)
{
  # run limma # 
  require(doParallel)
  require(foreach)
  require(iterators)    
  registerDoParallel(cores=cpus)    
  xlist <- split(1:reps, 1:cpus)
  foreach(i=1:cpus, .combine='cbind') %dopar% bootLimmaStat(FixedDataMatrix, CompleteTable, 
                                                            designTable, covarVec, 
                                                            j, xlist[[i]])
}


conf95 <- function(boot_table, i) 
{
  X <- mat.or.vec(nr=nrow(boot_table), nc=3)
  for (i in 1:nrow(boot_table)) {    
    ss <- na.omit(boot_table[i,])
    X[i,] <- quantile(ss, c(0.025, 0.5, 0.975))  
  }
  X <- as.data.frame(X)
  colnames(X) <- c("Left", "Med", "Right") 
  X
}


bootDiff <- function(designTable, dataMat, covarVec, writingDir, topTableCol=3, reps, cpus, ...){
  require(limma)
  
  FixedDataMatrix<-dataMat[, colnames(dataMat) %in% rownames(designTable)]
  
  #create file name for data output
  DataType       = strsplit(rownames(FixedDataMatrix)[1], split=":")[[1]][3]
  targetPheno <- unlist(str_split(colnames(designTable)[ncol(designTable)], ":"))[[5]]  # last col of design
  targetPheno <- str_replace_all(string = targetPheno, pattern = "`", replacement = "")
  print(targetPheno)
  OutputFile     = paste0("/DE_", DataType, "_", targetPheno)
  OutputFileName = add_date_tag(OutputFile, ".txt")
  OutputFile_dir = paste0(writingDir, OutputFileName)
  
  # make sure that the samples are the same in both the dataMatrix and design matrix
  CheckData<-table(colnames(FixedDataMatrix) %in% rownames(designTable))
  
  if (length(CheckData) == 1 & names(CheckData)[1] == "TRUE"){  
    Complete_table <- runLimma(FixedDataMatrix, designTable, covarVec)
    Boot_table <- runBootstrap(FixedDataMatrix, Complete_table, designTable, covarVec, cpus, topTableCol, reps)
    confInts <- conf95(Boot_table)
    SuperComplete <- cbind(Complete_table, confInts)
    write.table(SuperComplete, file=OutputFile_dir, quote = F, row.names = T)
    return(SuperComplete)
  }
  
  else{
    stop("Sample IDs in Design Matrix do not match Sample IDs in Data Matrix")
  }
}



add_date_tag <- function(stringToTag, fileExtension){
    today <- Sys.Date()
    todayf <- format(today, format="%Y%m%d")
    return(paste(stringToTag, "_", todayf, sep = "", fileExtension))
}
      
    
