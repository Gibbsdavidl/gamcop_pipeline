# General function to draw set size proporational pretty Venn Diagrams to the file given by the filename argument.
# need to explore parameters further -- right now can plot multicolored size proportional circles and save file 
# to a specific location (currently set all to defaults)
# max number of sets: 5 
# Author: Nyasha Chambwe
# Date: 01/28/2015

# arguments: 
# listOfItems: A list of vectors (e.g., integers, chars), with each component corresponding to a separate circle in the Venn diagram
# give each list item a name -> this will be the label for a given circle
# filename: Filename for image output, or if NULL returns the grid object itself
draw_venn_diagram <- function(listOfItems, filename){
  require(VennDiagram) || stop("Could not load package 'VennDiagram'")
  stopifnot(is.list(listOfItems))
  
  colors <- c( "seagreen3", "orchid3", "dodgerblue", "darkorange1", "green")
  
  venn.diagram(listOfItems, filename=filename, height = 3000, width = 4000, resolution = 500, 
               imagetype = "tiff", units = "px", compression = "lzw", na = "remove", main = "", sub = "", 
               main.pos = c(0.5, 1.05), main.fontface = "plain", main.fontfamily = "serif", 
               main.col = "black", main.cex = 1, main.just = c(0.5, 1), sub.pos = c(0.5, 1.05), 
               sub.fontface = "plain", sub.fontfamily = "serif", sub.col = "black", sub.cex = 1, 
               sub.just = c(0.5, 1), category.names = names(listOfItems), force.unique = TRUE, 
               fill = colors[1:length(listOfItems)]
  )
}

