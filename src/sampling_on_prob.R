find_probabilities<-function(target, orig_clin_Data, new_clin_Data){
  data_type<-strsplit(target, ":")
  data_type<-data_type[[1]][1]
  temp<-t(new_clin_Data[target,])
  items<-dat_stats(temp, data_type)

  categories<-dim(items)[1]
  
  alldat<-t(orig_clin_Data[target,])
  sample_probability<-alldat[,1]

  for (row in 2:categories){
    pheno<-items[row, 1]
    prob <-items[row, 3]
    sample_probability[ sample_probability==pheno ] <- prob
  }
  return(sample_probability)
}


dat_stats<-function(dat, data_type){
  vector<-data.frame(ID=character(), Samples=numeric(0), Prob=numeric(0), stringsAsFactors=FALSE)
  total_NA<-sum(is.na(dat))
  vector[1,]<-c("NA", total_NA, 0.00)
  
  if (data_type=="N"){
      vector[2,]<-c("NUM", sum(!is.na(dat)), 1.00)
  }else{
    a<-data.frame(table(dat))   
    total_samples<-sum(table(dat))
    for(row in 1:nrow(a)){
      a[,1]<-as.character(a[,1])
      row.no<-row+1
      prob<-a[row,2]/total_samples
      vector[row.no,]<-c(a[row,1], a[row,2], prob)
    }
  }
  return(vector)
}