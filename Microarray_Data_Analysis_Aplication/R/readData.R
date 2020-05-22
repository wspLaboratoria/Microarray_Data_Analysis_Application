readData <- function(dataPath, phenoDataPath){
  
  library(affy);library(gcrma); library(limma);library(oligo)
  
  annotatedDataFrame <- read.AnnotatedDataFrame(file.path(phenoDataPath),
                                                sep = "\t",
                                                header=TRUE,
                                                row.names = 4,
                                                stringsAsFactors = F)
  
  
  data <- read.affybatch(filenames = file.path(dataPath),
                         phenoData = annotatedDataFrame)
}