differential <- function(design, dataMat, covarVec, writingDir, ...){
   require(limma)
  
  FixedDataMatrix<-dataMat[, colnames(dataMat) %in% rownames(design)]
  
  #create file name for data output
  DataType       = strsplit(rownames(FixedDataMatrix)[1], split=":")[[1]][3]
  targetPheno <- unlist(str_split(colnames(design)[ncol(design)], ":"))[[5]]  # last col of design
  targetPheno <- str_replace_all(string = targetPheno, pattern = "`", replacement = "")
  print(targetPheno)
  OutputFile     = paste0("/DE_", DataType, "_", targetPheno)
  OutputFileName = add_date_tag(OutputFile, ".txt")
  OutputFile_dir = paste0(writingDir, OutputFileName)
  
  # make sure that the samples are the same in both the dataMatrix and design matrix
  CheckData<-table(colnames(FixedDataMatrix) %in% rownames(design))
  
  if (length(CheckData) == 1 & names(CheckData)[1] == "TRUE"){
    
    # run limma
    fit1<-lmFit(FixedDataMatrix, design)
    fit1<-eBayes(fit1)
    
    #gather relevant covariates
    cov_length = length(covarVec)+2
    DesignVariables  = dim(design)[2]
    Coef_of_Interest = paste(cov_length,DesignVariables, sep=":")
    
    #Write and Return Table of DE/M variables and the corresponding statistics
    print(cov_length:DesignVariables)
    Complete_table<-topTable(fit1, n=Inf, coef= cov_length:DesignVariables)
    write.table(Complete_table, file=OutputFile_dir, quote = F, row.names = T)
    return(Complete_table)
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
      
    
