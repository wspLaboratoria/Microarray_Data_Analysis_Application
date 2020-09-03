pca_genes <- function(dataNormalized){
  
  dataExprs <- as.data.frame(exprs(dataNormalized))
  colnames(dataExprs) <- dataNormalized@phenoData@data$CLASS
  pca <- prcomp(t(dataExprs))
  pc1 <-pca$x[,1]
  pc2 <-pca$x[,2]
  Klasy <- colnames(dataExprs)
  nazwy <- unique(dataNormalized@phenoData@data$CLASS)
  numeracja <- c(1:ncol(dataExprs))
  
  for(j in 1:length(nazwy)){
    for(i in 1:ncol(dataExprs)){
      if(Klasy[i]==nazwy[j]){
        numeracja[i] <- j}
    }
  }
  
  df_out <- as.data.frame(pca$x)
  p<-ggplot(df_out,aes(x=pc1,y=pc2,color=Klasy ))
  p<-p+geom_point()+ ggtitle("Wykres PCA")
  p
  return(p)
}