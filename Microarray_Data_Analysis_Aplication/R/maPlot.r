maPlot <- function(data){
  library(affy); library(oligo)
  
  if (class(data) == 'AffyBatch') {plot <- affy::MAplot(data,show.statistics = FALSE)}
  if (class(data) == 'ExpressionSet') {plot <- oligo::MAplot(data,plotFun=smoothScatter)}
  
  return(plot)
}