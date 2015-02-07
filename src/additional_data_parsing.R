##Parsing Data Scripts
##Created by: Crystal Humphries
##Objective: to streamline data parsing based on clinical data 

# returns colnames of 'X.101.443.FAM' to 101-443-FAM
# e.g. colnames(df)<-cor_colnames(df)

cor_colnames<-function( dat.Frame){
  col.dat<-colnames(dat.Frame)
  col.dat<-gsub("X", "", col.dat)
  col.dat<-gsub("\\.", "-", col.dat)
  return(col.dat)
}



#subsets data based on ONE trait

data_subset_Clin<-function(Mol.dat.frame, Clin.dat.frame, clin.rowname, traits){
    #clin.rowname = rownames(Clin.dat)['desired trait']
    #traits        = vector of matching items wanted (e.g. blood draw date ranges from -31  to 309 days. To get days 1 and 2, 
    # make a vector of traits = c(1,2) )
 
  nm  = row.names(Clin.dat) %in% clin.rowname
  if (sum(nm) ==0){
    print("data_subset_Clin error: No Matching Clinical Trait")
    return(NA)
  }
  
  desired_trait  = t(Clin.dat[ nm,])
  desired_trait = data.frame(desired_trait)
  
  names(desired_trait) = "Trait"
  if (substr(clin.rowname,1,1) == "N"){
    desired_trait$Trait <-as.numeric(as.character(desired_trait$Trait))
  }
  
  wanted_samples<-(desired_trait$Trait %in% traits)
  
  if (sum(wanted_samples) == 0){
    print("data_subset_Clin error: No samples to subset")
    return(NA)
  }
  
  desired_trait<-subset.data.frame(desired_trait, wanted_samples)
  if (dim(desired_trait)[1] == 0){
    print("data_subset_Clin error: No data to subset")
    return(NA)
  }
    dat.frame.new<-(dat.frame[,(colnames(dat.frame) %in% rownames(desired_trait))])
    return(dat.frame.new)
}