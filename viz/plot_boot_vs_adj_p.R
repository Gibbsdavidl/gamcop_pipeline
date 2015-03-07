

# plotting the adjusted p-values vs the bootstrap interval

plotBootVsAdjP <- function(deTable, atitle){  
  plot(x=c(0,1), y=c(0,1), pch='$', xlab = "Adjusted P-Value", ylab="Bootstrap adj. p-value")
  abline(b=1,a=0)
  points(deTable$adj.P.Val, Preterm$Left, col="red")
  points(deTable$adj.P.Val, Preterm$Right, col="blue")
  points(deTable$adj.P.Val, Preterm$Med)
}
