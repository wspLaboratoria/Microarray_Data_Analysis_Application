# normalization <- function(microArrayData, phenoData, normType) {
normalization <- function(data, normType) {
  
  library(affy);library(gcrma); library(limma);library(oligo)
  
  # annotatedDataFrame <- read.AnnotatedDataFrame(file.path(phenoData),
  #                                               sep = "\t",
  #                                               header=TRUE,
  #                                               row.names = 4,
  #                                               stringsAsFactors = F)
  # 
  # 
  # data <- read.affybatch(filenames = file.path(microArrayData),
  #                                     phenoData = annotatedDataFrame)
  
  if(normType == "RMA"){normData = affy::rma(data)}
  if(normType == "GCRMA"){normData = gcrma(data)}
  if(normType == "MAS"){normData = mas5(data)}

  return(normData)

}
