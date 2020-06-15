histPlot=function(data,data.norm,geneNumber){
  
  data.matrix=exprs(data)
  data.matrix1=exprs(data.norm)
  
  t <- list(size = 11)
  
  if(identical(data.matrix,data.matrix1)==T){
    
    
    wariancje <- as.matrix(apply(data.matrix1,1,var))
    wariancjeIndeksyPosortowane <- order(wariancje,decreasing = TRUE)
    
    data.matrix1 <- data.matrix1[wariancjeIndeksyPosortowane[1:geneNumber],]
    
    sampleNames =  vector() 
    logs =  vector() 
    for  (i in  1 : dim(data.matrix)[ 2 ]) { 
      sampleNames =  c(sampleNames, rep(colnames(data.matrix)[i] , dim(data.matrix)[1])) 
      logs =  c(logs,(data.matrix[,i]))
    }
    
    logData =  data.frame(logInt=logs,sampleName=sampleNames )
    
    dataHist2 = ggplot(logData,aes(logInt,color=sampleName))+geom_density ()

    fig=ggplotly(dataHist2)
    fig=layout(fig, title=c('Histogram'),yaxis=list(title='density'),
               xaxis=list(title='logInt'),font=t)
    fig
    
  }else{
    
    wariancje <- as.matrix(apply(data.matrix1,1,var))
    wariancjeIndeksyPosortowane <- order(wariancje,decreasing = TRUE)
    
    data.matrix1<- data.matrix1[wariancjeIndeksyPosortowane[1:geneNumber],]
    genes=rownames(data.matrix1)
    
    
    
    tmp=c()
    for(i in genes){
      tmp=rbind(tmp,colMeans(affy::pm(data,i)))
    }
    
    
    sampleNames =  vector () 
    logs =  vector () 
    for  (i in  1 : dim(tmp)[2]) { 
      sampleNames =  c (sampleNames, rep(colnames(tmp)[i], dim (tmp)[1])) 
      logs =  c (logs, log(tmp[,i]))
    }
    
    
    logData =  data.frame(logInt = logs, sampleName = sampleNames)
 
    
    
    dataHist2 = ggplot(logData, aes(logInt,color=sampleName))+ geom_density ()
    
    fig=ggplotly(dataHist2)
    fig=layout(fig, title=c('Histogram'),yaxis=list(title='density'),
               xaxis=list(title='logInt'),font=t)
    fig

    
  }
  
  return(fig)
  
}