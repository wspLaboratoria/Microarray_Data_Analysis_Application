dataTableDifferencialGenes <- function(dataNormalized,grupa1,grupa2){
  dataExprs <- exprs(dataNormalized)
  colnames(dataExprs) <- dataNormalized@phenoData@data$CLASS
  kolumny <- colnames(dataExprs)
  data1 <- dataExprs[,kolumny==grupa1]
  data2 <- dataExprs[,kolumny==grupa2]
  geny <-rownames(dataExprs)
  results <- c()
  for(i in 1:nrow(dataExprs)){
    tmp <- wilcox.test(data1[i,],data2[i,],alternative = "two.sided")
    tmp <- tmp$p.value
    wiersz <- c(geny[i],tmp)
    results <- rbind(results,wiersz)
  }
  as.data.frame(results)
  colnames(results)<-c("Gen","p_value")
  rownames(results) <- c()
  
  #poprawka FDR
  FDR_p_values <- as.matrix(p.adjust(results[,2],method='fdr',n=length(results[,2])))
  colnames(FDR_p_values) <- "p_value FDR"
  
  results <- cbind(results,FDR_p_values)
  
  #roznicujace
  Table_Diff_Genes <- as.data.frame(results[results[,"p_value FDR"]<0.05,])
  return(Table_Diff_Genes)
}