boxplot=function(data,data.norm,geneNumber){
  
  data.matrix=exprs(data)
  data.matrix1=exprs(data.norm)
  
  t <- list(size = 8)
  
  if(identical(data.matrix,data.matrix1)==T){
    
    # colnames(data.matrix1)=description@data[["Sample"]]
    
    wariancje <- as.matrix(apply(data.matrix1,1,var))
    wariancjeIndeksyPosortowane <- order(wariancje,decreasing = TRUE)
    
    data.matrix<- data.matrix1[wariancjeIndeksyPosortowane[1:geneNumber],]
    
  
    fig <- plot_ly(x = data.matrix1[,1], type = "box",name=colnames(data.matrix1)[1])
    
    for(i in 2:dim(data.matrix1)[2]){
      
      fig <- fig %>% add_trace(x = data.matrix1[,i],name=colnames(data.matrix1)[i])
    }
    fig=layout(fig, title=c('Boxplot'),yaxis=list(title='Sample Names',zeroline=F,showline=T,mirror='ticks'),
               xaxis=list(title='logInt',zeroline=T,showline=T,mirror='ticks',font=t))
    
    fig
    
    
  }else{
    
    # colnames(data.matrix1)=description@data[["Sample"]]
    
    wariancje <- as.matrix(apply(data.matrix1,1,var))
    wariancjeIndeksyPosortowane <- order(wariancje,decreasing = TRUE)
    
    data.matrix1 <- data.matrix1[wariancjeIndeksyPosortowane[1:geneNumber],]
    genes=rownames(data.matrix1)
  
    
    tmp=c()
    for(i in genes){
      tmp=rbind(tmp,colMeans(affy::pm(data,i)))
    }
   
    tmp=log(tmp)
    
    fig <- plot_ly(x = tmp[,1], type = "box",name=colnames(tmp)[1])
    
    for(i in 2:dim(tmp)[2]){
      
      fig <- fig %>% add_trace(x = tmp[,i],name=colnames(tmp)[i])
    }
    fig=layout(fig, title=c('Boxplot'),yaxis=list(title='Sample Names',zeroline=F,showline=T,mirror='ticks'),
               xaxis=list(title='logInt',zeroline=T,showline=T,mirror='ticks',font=t))
    
    fig

  }
  
  return(fig)
  
}
