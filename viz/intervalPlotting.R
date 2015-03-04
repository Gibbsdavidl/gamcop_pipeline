

# plotting some intervals.

Preterm <- read.delim("/Volumes/StorageDisk/Meth_DF5/pipeline/DEs_Feb_27th/DE_blood_day_04/Preterm.txt", stringsAsFactors=FALSE)
Idiop <- read.delim("/Volumes/StorageDisk/Meth_DF5/pipeline/DEs_Feb_27th/DE_blood_day_04/Inova_Idiopathic_NA.txt", stringsAsFactors=FALSE)


drawInterval <- function(deTable, rowI, atitle) {
  
  adjP  <- deTable$adj.P.Val[rowI]
  bootL <- deTable$Left[rowI]
  bootM <- deTable$Med[rowI]
  bootR <- deTable$Right[rowI]
  
  seg1   <- min(c(adjP, bootL, bootM, bootR))
  seg2   <- max(c(adjP, bootL, bootM, bootR))
  segby  <- abs(seg2 - seg1) / 10 
  range1 <- min(c(adjP, bootL, bootM, bootR)) - 5*segby
  range2 <- max(c(adjP, bootL, bootM, bootR)) + 5*segby
  
  plot.new() 
  plot.window( c(range1, range2), c(-1, 1), main=atitle) 
  title(atitle)

  ticksX <- signif(seq(from = range1 , to= range2, by=segby), 2)
  
  axis(1, at=ticksX, pos=0) 

  # adj. p-value
  lines( c(bootL,bootL,bootR,bootR), c( -0.25, -0.5, -0.5, -0.25 ) ) 
  text( bootL+segby, -0.6, '95% Bootstrap Interval' )  # the plotrix package has a function for text in a box 

  lines( c(adjP,adjP), c( 0.1, 0.3) ) 
  text( adjP, 0.4, 'Adj. P-Value' ) 

  lines( c(0.05,0.05), c( 0.1, 0.5) ) 
  text( 0.05, 0.55, '5% Signif. Threshold' ) 

  if (range2 > 0.2) {
    lines( c(0.2,0.2), c( 0.1, 0.75) ) 
    text( 0.2, 0.75, '20% Signif. Threshold' ) 
  }
  
  lines( c(bootM,bootM), c( -0.1, -0.8) ) 
  text( bootM, -0.9, 'Median Bootstrap Adj. P-value' ) 
  
  points(x = 0, y=0, lwd=2, col="red")
  
}



