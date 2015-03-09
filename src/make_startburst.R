#3/2/2015
#by Crystal Humphries


#plot for visualizing Overlapping data in a starburst plot
#include a list of p-values and logFC. 
#uses output generated from the "Gene_level_overlap.R" function
#generate list in the following manner: list(RNA=c(logFC_cutoff, Pvalue_cutoff), Meth=c(logFC_cutoff, Pvalue_cutoff))
##############
#
# the default FDR cutoff is 0.05 
# the default logFC cutoff is 0.5
####################

require(ggplot2)
require(reshape2)
require(RColorBrewer)

make_starburst<-function(dataset1, dataset2, infoList, writingDir){
  #check NameVector
  NameVector = names(infoList)
  if(length(NameVector)!=2){
      print ("infoList is not complete")
      return(NA)
  }
  if (is.na(data.list[[1]][2]) | !(is.numeric(data.list[[1]][2]))){
    data.list[[1]][2]<-0.05
  }
  if (is.na(data.list[[2]][2]) | !(is.numeric(data.list[[2]][2]))){
    data.list[[2]][2]<-0.05
  }
  if (is.na(data.list[[1]][1]) | !(is.numeric(data.list[[1]][1]))){
    data.list[[1]][1]<-0.5
  }
  if (is.na(data.list[[2]][1]) | !(is.numeric(data.list[[2]][1]))){
    data.list[[2]][1]<-0.5
  }
  
  merged_table<-generate_de_overlap_table(dataset1, dataset2)
  names(merged_table)<-gsub("x", NameVector[1], names(merged_table))
  names(merged_table)<-gsub("y", NameVector[2], names(merged_table))
  
  dataset<-NULL
  for(i in 1:2){
    dataset<-c(dataset,paste("logFC",NameVector[i], sep="."))
    dataset<-c(dataset,paste("adj.P.Val",NameVector[i], sep="."))
  }
  a<-which(names(merged_table) %in% dataset)
  
  if(length(a)!=4){
      print("Colnames are not from Limma Standard output. Make sure each data set has a adj.P.Val and a logFC column")
      return("NA")
  }
  
  starburst.m<-merged_table[,a]

  #color differentially expressed genes as dots
  PvalThres.1<-infoList[[1]][2]
  PvalThres.2<-infoList[[2]][2]
  
  #get FC cutoff lines ( be sure to include both negeative and positive values)
  FCcutOffs.1<-c(infoList[[1]][1], -infoList[[1]][1])
  FCcutOffs.2<-c(infoList[[2]][1], -infoList[[2]][1]) 
  
  SuperSig.x<-(abs(starburst.m[,1])>FCcutOffs.1[1] & starburst.m[,2]<PvalThres.1)
  SuperSig.y<-(abs(starburst.m[,3])>FCcutOffs.2[1] & starburst.m[,4]<PvalThres.2)
  starburst.m<-data.frame(starburst.m, SuperSig.x, SuperSig.y) 

  Sig<-c(rep("Not Significant", dim(starburst.m)[1]))
  Sig[ starburst.m$SuperSig.x=="TRUE" & starburst.m$SuperSig.y=="FALSE"]<-paste0("Signif_",NameVector[1])
  Sig[ starburst.m$SuperSig.x=="FALSE" & starburst.m$SuperSig.y=="TRUE"]<-paste0("Signif_",NameVector[2])
  Sig[ starburst.m$SuperSig.x=="TRUE" & starburst.m$SuperSig.y=="TRUE"]<-paste0("Signif_",NameVector[1], "&", NameVector[2])
  starburst.m<-data.frame(starburst.m, Sig)
  
  names(starburst.m)<-gsub("x", NameVector[1], names(starburst.m))
  names(starburst.m)<-gsub("y", NameVector[2], names(starburst.m))
  
 cbPalette<-c("#404040",  "#377eb8","#4daf4a", "#e41a1c")
    
 p<- ggplot(starburst.m, aes(x=starburst.m[,1], y=starburst.m[,3], colour=Sig)) + geom_point()  
 p<- p + geom_vline(xintercept=0, size=0.5) + geom_hline(yintercept=0, size=0.5)  
 p<- p + geom_hline(yintercept=FCcutOffs.2[1], linetype="dotted") + geom_hline(yintercept=FCcutOffs.2[2], linetype="dotted") 
 p<- p + xlab(paste(NameVector[1], "\n", paste0("|LogFC|>", FCcutOffs.1[1], " and FDR<", PvalThres.1))) + ylab(paste(NameVector[2], "\n", paste0("|LogFC|>", FCcutOffs.1[2], " and FDR<", PvalThres.2)))
 p<- p + ggtitle(paste("Overlap of", NameVector[1], "and", NameVector[2])) + scale_colour_manual(values=cbPalette) 
 p<- p + geom_vline(xintercept=FCcutOffs.1[1], linetype="dotted") + geom_vline(xintercept=FCcutOffs.1[2], linetype="dotted") 
 p<- p + theme(plot.title=element_text(face="bold", size=20), legend.title=element_blank(), legend.position="top", legend.text= element_text(size=14)) 
 p<- p + theme(axis.title.y=element_text(size=16, fac="bold"), axis.title.x=element_text(size=16, face="bold"))
 print(p)
 
 file_name.beta<-paste0("Starburst.", NameVector[1], "_",NameVector[2])
 file_name<-add_date_tag(file_name.beta, ".pdf")
 file<-paste0(writingDir, file_name)
 ggsave(p, file=file)
}

add_date_tag <- function(stringToTag, fileExtension){
  today <- Sys.Date()
  todayf <- format(today, format="%Y%m%d")
  return(paste(stringToTag, "_", todayf, sep = "", fileExtension))
}
