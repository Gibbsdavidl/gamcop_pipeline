# The set of plotting functions to visualize differential expression from Limma differential 
# expression analysis tables

# General Helper Functions ------------------------------------------------
# Function to add a date tag to the end of a file name
# @param stringToTag string to append date to
# @param fileExtension to put at the end
add_date_tag <- function(stringToTag, fileExtension){
  today <- Sys.Date()
  todayf <- format(today, format="%Y%m%d")
  return(paste(stringToTag, "_",todayf, sep = "", fileExtension))
}
# Function to beautify ggplots in a generic way - changing background to black and white, axis labels and main titile for each figure
makeNeatGraphs <- function(plotObj, xlab, ylab, newTitle){
  
  # If need to change some of the colors
  # The palette with grey:
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  # The palette with black:
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  
  plotObj <- plotObj + 
    labs(title=newTitle, x= xlab,  y=ylab) +
    theme_bw() + scale_fill_manual(values=cbPalette) +
    theme(axis.title.x= element_text(face= "bold", size=20, colour="black", vjust = 0.3)) + 
    theme(axis.title.y= element_text(face= "bold", size=20, colour="black", hjust = 0.3)) +
    theme(axis.text.x=element_text(face="bold", size=16, colour="black", vjust = 0.5)) +
    theme(axis.text.y=element_text(face="bold", size=16, colour="black")) +
    theme(legend.justification=c(1,1), legend.position ="top") +
    theme(title=element_text(size=20, face="bold"))
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
  
  filename <- add_date_tag(paste(writingDir, "/volcano_plot", comparisonName, dataSource, sep = "_"), fileExtension = ".pdf")
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
expression_by_phenotype_boxplots <- function(dat, gene, pheno){
  print(expression_by_phenotype_boxplots(dat, gene_of_interest, phenotype_of_interest))
  names(dat) <- c("Phenotype_class", "miRNA_expression")
  plot.obj <- ggplot(dat, aes(x = Phenotype_class, y=miRNA_expression, fill=Phenotype_class) )
  plot.obj <- plot.obj + geom_boxplot() + labs(title=gene) + theme(axis.title.x = element_blank(), legend.position="top") +   # Remove x-axis label
    ylab("Log2 Normalized Expression")  
  plot.obj <- makeNeatGraphs(plotObj = plot.obj, xlab = "", ylab="Log2 Normalized Expression", newTitle = gene)
  
  filename <- add_date_tag(paste("results/box-plots-", comparisonName, gene_of_interest, sep=""),fileExtension = ".pdf")
  pdf(file = filename , width = 8, height = 6)
  print(plotObj)
  dev.off()
}

# Heatmap Part ------------------------------------------------------------
