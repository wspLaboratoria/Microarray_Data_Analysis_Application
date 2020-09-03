readData <- function(dataPath, phenoDataPath){
  
  
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
