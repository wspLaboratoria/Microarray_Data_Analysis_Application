boxplot = function (data){
  
  data.matrix=exprs(data)
  fig <- plot_ly(x = data.matrix[,1], type = "box",name=colnames(data.matrix)[1])
  
  for(i in 2:dim(data.matrix)[2]){
    
    fig <- fig %>% add_trace(x = data.matrix[,i],name=colnames(data.matrix)[i])
  }
  fig=layout(fig, title=c('Tytu? '),yaxis=list(title='Sample Names',zeroline=F,showline=T,mirror='ticks'),
             xaxis=list(title='O? X',zeroline=T,showline=T,mirror='ticks',font=t))
  
  fig
  
}