# The set of plotting functions to visualize differential expression from Limma differential 
# expression analysis tables

# General Helper Functions ------------------------------------------------
# Function to add a date tag to the end of a file name
# @param stringToTag string to append date to
# @param fileExtension to put at the end
add_date_tag <- function(stringToTag, fileExtension){
  today <- Sys.Date()
  todayf <- format(today, format="%Y%m%d")
  return(paste(stringToTag, "_", todayf, sep = "", fileExtension))
}
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

# Function to beautify ggplots in a generic way - changing background to black and white, axis labels and main titile for each figure
makeNeatGraphs <- function(plotObj, xlab, ylab, newTitle, legendTitle){
  
  # If need to change some of the colors
  # The palette with grey:
  cbPalette <- c( "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#999999","#D55E00", "#CC79A7",  "#F0E442")
  
  # The palette with black:
  cbbPalette <- c( "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#000000", "#D55E00", "#CC79A7",  "#F0E442")
  
  
  plotObj <- plotObj + 
    labs(title=newTitle, x= xlab,  y=ylab) +
    theme_bw() + scale_fill_manual(name=legendTitle, values=cbPalette) +
    theme(axis.title.x= element_text(face= "bold", size=16, colour="black", vjust = 0.3)) + 
    theme(axis.title.y= element_text(face= "bold", size=16, colour="black", hjust = 0.3)) +
    theme(axis.text.x=element_text(face="bold", size=14, colour="black", vjust = 0.5)) +
    theme(axis.text.y=element_text(face="bold", size=14, colour="black")) +
    theme(legend.justification=c(1,1), legend.position ="top") +
    theme(title=element_text(size=18, face="bold"))
  return(plotObj)
}

# Volcano Plot ------------------------------------------------------------
# Takes the output of Limma topTable and returns a volcano plot 
# with -log10(Adj.P.Value) on the y-axis and log2(FC) on the x-axis
# points are colored:
#   red: if they pass the 
#   orange: if they pass the fold-change threshold 
# output of limma is augmented with a column "Gene" with labels for points
make_volcano_plot <- function(topTable, dataSource, pValueThresh, fcThreshold, comparisonName, writingDir){
  # code source: https://gist.github.com/stephenturner/4a599dbf120f380d38e7#file-volcanoplot-r
  require(calibrate) || stop("Could not load package 'calibrate'")
  foldChangeLim <-  round(max(abs(range(topTable$logFC)))) + 0.5
  
  ouputFileName <- paste0(writingDir,"/","volcano_plot_", comparisonName, dataSource, sep = "_")
  filename <- add_date_tag(ouputFileName, fileExtension = ".pdf")
  pdf(file = filename , width = 8, height = 6)
  
  # Make a basic volcano plot
  with(topTable,{ 
    plot(logFC, -log10(adj.P.Val), pch=20, col="gray", bg="gray", 
         main=comparisonName,
         xlab="log2(Fold Change)", ylab="-log10(B.H. Adjusted P value)", xlim=c(-foldChangeLim, foldChangeLim))
    abline(v=fcThreshold, col="gray", lty=3)
    abline(v=-(fcThreshold),col="gray",lty=3)
    abline(h=-log10(pValueThresh), col="gray", lty=3)
  })
  
  # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
  with(subset(topTable, adj.P.Val< pValueThresh ), points(logFC, -log10(adj.P.Val), pch=20, col="red"))
  with(subset(topTable, abs(logFC)> fcThreshold), points(logFC, -log10(adj.P.Val), pch=20, col="orange"))
  with(subset(topTable, adj.P.Val<pValueThresh & abs(logFC)>fcThreshold), points(logFC, -log10(adj.P.Val), pch=20, col="green"))
  
  # Label points with the textxy function from the calibrate plot
  # with(subset(topTable, adj.P.Val<pValueThresh & abs(logFC)>fcThreshold), textxy(logFC, -log10(adj.P.Val), labs=Gene, cex=.8, offset=0.8))
  dev.off()
  
}

# Gene level box plots ----------------------------------------------------
# Boxplots for a given gene
# phenotype must be a factor
expression_by_phenotype_boxplots <- function(dat, gene, pheno, dataSource, writingDir){
  names(dat) <- c("PhenotypeClass", "Expression")
  dat$PhenotypeClass <- as.factor(dat$PhenotypeClass)
  plot.obj <- ggplot(dat, aes(x = PhenotypeClass, y=Expression, fill=PhenotypeClass) )
  plot.obj <- plot.obj + geom_boxplot(width=0.5) 
  if(dataSource=="RNASeq" || dataSource=="MIRNA"){
    plot.obj <- makeNeatGraphs(plotObj = plot.obj, xlab = "", ylab="Log2 Normalized Expression", newTitle = gene, pheno)
  } else{
    plot.obj <- makeNeatGraphs(plotObj = plot.obj, xlab = "", ylab="Log2 Normalized Methylation", newTitle = gene, pheno)
  }
  plot.obj <- plot.obj + labs(title=gene) + theme(axis.title.x = element_blank(), legend.position="top")  # Remove x-axis label
  filename <- add_date_tag(paste(paste0(writingDir,"/" ,"box_plot"), dataSource, pheno, gene, sep="_"),fileExtension = ".pdf")
  pdf(file = filename , width = 8, height = 6)
  print(plot.obj)
  dev.off()
}

# Heatmap Part ------------------------------------------------------------