readData <- function(dataPath, phenoDataPath){
  
  library(affy);library(gcrma); library(limma);library(oligo)
  
  annotatedDataFrame <- read.AnnotatedDataFrame(file.path(phenoDataPath),
                                                sep = "\t",
                                                header=TRUE,
                                                row.names = 4,
                                                stringsAsFactors = F)
  
  
  data <- ReadAffy(filenames = file.path(dataPath),
                   sampleNames = annotatedDataFrame@data[["Sample"]],
                   phenoData = annotatedDataFrame)


  return(data)
}
