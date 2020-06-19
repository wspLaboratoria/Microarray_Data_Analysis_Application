kmen <- function(data, geneNumber, number_klast, dist_method){
  
  
  # library(tidyverse)  # data manipulation
  # library(cluster)    # clustering algorithms
  # library(factoextra)
  # library(amap)

  data.matrix = exprs(data)
  wariancje <- as.matrix(apply(data.matrix,1,var))
  wariancjeIndeksyPosortowane <- order(wariancje,decreasing = TRUE)
  
  dataHeatmap <- data.matrix[wariancjeIndeksyPosortowane[1:geneNumber],]
  mydatascale <- (scale(t(dataHeatmap))) # Centers
  
  
  k2 <- Kmeans(mydatascale, centers = number_klast, method = dist_method)
  
  
  plot_km <- fviz_cluster(k2, data = mydatascale, geom='text')
  
  tmp=data.frame(Clusters=k2$cluster)
  table=data.frame(Sample=rownames(tmp),tmp,row.names = NULL) 
    
    
  
  return(list(plot_km,table))

  
}