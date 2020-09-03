kmen=function(data,geneNumber,number_klast,dist_method){
  
  
  # library(tidyverse)  # data manipulation
  # library(cluster)    # clustering algorithms
  # library(factoextra) 
  # library(amap)
  # library(colorspace)
  
  data.matrix = exprs(data)
  wariancje <- as.matrix(apply(data.matrix,1,var))
  wariancjeIndeksyPosortowane <- order(wariancje,decreasing = TRUE)
  
  dataHeatmap <- data.matrix[wariancjeIndeksyPosortowane[1:geneNumber],]
  mydatascale <- (scale(t(dataHeatmap))) # Centers
  
  
  k2 <- Kmeans(mydatascale, centers = number_klast, method = dist_method)
  rc <- colorspace::rainbow_hcl(number_klast)
  
  
  
  plot_km=fviz_cluster(k2, data = mydatascale, palette = rc)
  
 
  p4=ggplotly(plot_km)
  
  for (i in 2:(length(p4[["x"]][["data"]])/2) ){
    p4[["x"]][["data"]][[i]][["name"]] <-str_remove_all(p4[["x"]][["data"]][[i]][["name"]],'[(,1,NA)]')
    p4[["x"]][["data"]][[i]][["legendgroup"]]=str_remove_all(p4[["x"]][["data"]][[i]][["name"]],'[(,1,NA)]')
    
  }
  for (i in ((length(p4[["x"]][["data"]])/2)+1): length(p4[["x"]][["data"]])){
    p4[["x"]][["data"]][[i]][["showlegend"]]=F
    p4[["x"]][["data"]][[i]][["legendgroup"]]=str_remove_all(p4[["x"]][["data"]][[i]][["name"]],'[(,1,NA)]')
    p4[["x"]][["data"]][[16]][["legendgroup"]]='1'; p4[["x"]][["data"]][[16]][["name"]]='1'
    #   colors=c()
    #for( j in 1:length(p4[["x"]][["data"]][[i]][["text"]])){
    
    #   if ( str_sub(p4[["x"]][["data"]][[i]][["text"]][j],1,1)==str_sub(unique(data@phenoData@data$CLASS)[1],1,1)){
    
    #    colors=c(colors,'red')
    # }else if( str_sub(p4[["x"]][["data"]][[i]][["text"]][j],1,1)==str_sub(unique(data@phenoData@data$CLASS)[2],1,1)){
    #     colors=c(colors,'green')
    #  }else if( str_sub(p4[["x"]][["data"]][[i]][["text"]][j],1,1)==str_sub(unique(data@phenoData@data$CLASS)[3],1,1)){
    #   colors=c(colors,'blue')
    #}else if( str_sub(p4[["x"]][["data"]][[i]][["text"]][j],1,1)==str_sub(unique(data@phenoData@data$CLASS)[4],1,1)){
    #  colors=c(colors,'black')
    #}else if( str_sub(p4[["x"]][["data"]][[i]][["text"]][j],1,1)==str_sub(unique(data@phenoData@data$CLASS)[5],1,1)){
    # colors=c(colors,'magenta')
    #}
    
    #}
    #  p4[["x"]][["data"]][[i]][["textfont"]][["color"]]=colors
    
  }
  
  plot_km=p4
  

  
tmp=data.frame(Clusters=k2$cluster)
table=data.frame(Sample=rownames(tmp),tmp,row.names = NULL) 
  
  

return(list(plot_km,table))

  
  
}
