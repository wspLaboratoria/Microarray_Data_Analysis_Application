BoxPlotDifferencialGenes <- function(dataNormalized,grupa1,grupa2,progFDR){

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
  colnames(results)<-c("Gene","P value")
  rownames(results) <- c()
  
  #poprawka FDR
  FDR_p_values <- as.matrix(p.adjust(results[,2],method='fdr',n=length(results[,2])))
  colnames(FDR_p_values) <- "P value FDR"
  
  results <- cbind(results,FDR_p_values)
  
  #roznicujace
  Table_Diff_Genes <- as.data.frame(results[results[,"P value FDR"]<progFDR,])
  
  if(nrow(Table_Diff_Genes) ==0 ){ fig<- plot.new()}
  else{
  dataToPlot <- data.frame(0,1,3)
  for(i in 1:nrow(Table_Diff_Genes)){
    tmp_gene <- as.character(Table_Diff_Genes$Gene[i])
    tmp_1 <- as.data.frame(data1[tmp_gene,])
    tmp_2 <- as.data.frame(data2[tmp_gene,])
    tmp_1 <- cbind(tmp_1,replicate(nrow(tmp_1), grupa1),tmp_gene)
    tmp_2 <- cbind(tmp_2,replicate(nrow(tmp_2), grupa2),tmp_gene)
    colnames(tmp_1) <- colnames(tmp_2) <- colnames(dataToPlot) <- c("Ekspresja","Grupa","Gen")
    dataToPlot <- rbind(dataToPlot,tmp_1,tmp_2)
    }
    
  dataToPlot <- dataToPlot[2:nrow(dataToPlot),]
  rownames(dataToPlot) <- c()
  
  fig <- plot_ly(dataToPlot, x = ~Gen, y = ~Ekspresja, color = ~Grupa, type = "box")
  fig <- fig %>% layout(boxmode = "group")
  fig <- layout(fig, title=c(paste("Box plot for differential genes - two groups. FDR: ",progFDR)))
  }
    
  return(fig)
}