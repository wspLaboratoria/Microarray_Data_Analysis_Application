
createBoxplot <- function(DATA){

  library(affy);library(gcrma); library(limma)  
  plot = boxplot(DATA)
  
  return(plot)
}
