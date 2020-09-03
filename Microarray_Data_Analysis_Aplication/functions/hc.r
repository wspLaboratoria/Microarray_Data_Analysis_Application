hc=function(data,geneNumber,dist_method,clust_method,cut_number){
  
  # library(dendextend)
  # library(factoextra)
  # library(ggdendro)
  # library(colorspace)
  # library(heatmaply)
  # library(stringr)
  #library(gplots) jak zwykla hitmapa 
 # cut_number=5
 # geneNumber=10
  #clust_method='complete'
  #dist_method='euclidean'
  
  
  
  #dist_method= "euclidean", "maximum", "manhattan", "canberra"
  #clust_method='average','complete','single','centroid'
  
  data.matrix = exprs(data)
  wariancje <- as.matrix(apply(data.matrix,1,var))
  wariancjeIndeksyPosortowane <- order(wariancje,decreasing = TRUE)
  
  dataHeatmap <- data.matrix[wariancjeIndeksyPosortowane[1:geneNumber],]
  #head(mydatascale)
  
  
  mydatascale <- (scale(t(dataHeatmap))) # Centers
  d <- dist(mydatascale, method = dist_method)
  hc <- hclust(d, method =clust_method)
  
  rc <- colorspace::rainbow_hcl(cut_number)
  
 # plot(hc)
  sub_grp <- cutree(hc, k = cut_number)
  
  klastry=fviz_cluster(list(data = mydatascale, cluster = sub_grp), palette = rc)
  
  
  
  klastry
  p4=ggplotly(klastry)
  
  # hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(length(unique(data$CLASS)))
  
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
    
  klastry=p4
  
  paleta=c()
  for(w in 1:cut_number){
    paleta=c(paleta,(rc[w]))
    
    
  }
  names(paleta)=c(1:cut_number)
 
  p <- heatmaply(dataHeatmap, 
                 dendrogram = c("both"),
                 dist_method = dist_method,
                 hclust_method = clust_method,
                 show_dendrogram = c(TRUE, TRUE),
                 Colv = hc,
                 row_dend_left=F,
                 label_names = c("Gene", "Sample",'value'),
                 margins = c(60,100,40,20),
                 #grid_color = "white",
                 #grid_width = 0.00001,
                 hide_colorbar = T,
                 branches_lwd = 0.1,
                 fontsize_row = 5, fontsize_col = 5,

                 col_side_colors = sub_grp,
                 col_side_palette = paleta,
                 trace='none'
  )


  return(list(klastry,p))

  
  
  
  
}
