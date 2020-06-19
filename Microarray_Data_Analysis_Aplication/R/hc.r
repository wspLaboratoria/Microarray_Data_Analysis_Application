hc=function(data,geneNumber,dist_method,clust_method,cut_number){
  
  # library(dendextend)
  # library(factoextra)
  # library(ggdendro)
  # library(colorspace)
  # library(heatmaply)
  
  # dist_method= "euclidean", "maximum", "manhattan", "canberra"
  # clust_method='average','complete','single','centroid'
  
  data.matrix = exprs(data)
  wariancje <- as.matrix(apply(data.matrix,1,var))
  wariancjeIndeksyPosortowane <- order(wariancje,decreasing = TRUE)
  
  dataHeatmap <- data.matrix[wariancjeIndeksyPosortowane[1:geneNumber],]
  #head(mydatascale)
  
  
  mydatascale <- (scale(t(dataHeatmap))) # Centers
  d <- dist(mydatascale, method = dist_method)
  hc <- hclust(d, method =clust_method)
  
  
  # plot(hc)
  sub_grp <- cutree(hc, k = cut_number)
  
  klastry <- fviz_cluster(list(data = mydatascale, cluster = sub_grp), geom='text')
  
  
  #heatmap.2(dataHeatmap, 
  #         dendrogram='col', 
  #Rowv =as.dendrogram(hr), Colv= as.dendrogram(hc),
  #        label_names = c("row", "column", "value"),
  #       col = viridis(100),
  #      trace="none",              # hide trace
  #     density.info="none",       # hide histogram
  
  #    margins = c(5,18),         # margin on top(bottom) and left(right) side.
  #   cexRow=1, cexCol = 0.8,      # size of row / column labels
  #  srtCol=0, adjCol = c(0.5,1), # adjust the direction of row label to be horizontal
  # margin for the color key
  # ("bottom.margin", "left.margin", "top.margin", "left.margin" )
  #key.par=list(mar=c(5,1,3,1))
  # key = F,
  # )
  p <- heatmaply(dataHeatmap, 
                 dendrogram = c("both"),
                 dist_method = dist_method,
                 hclust_method = clust_method,
                 show_dendrogram = c(TRUE, TRUE),
                 Colv = hc,
                 row_dend_left=F,
                 label_names = c("row", "column",'value'),
                 margins = c(60,100,40,20),
                 #grid_color = "white",
                 #grid_width = 0.00001,
                 #titleX = FALSE,
                 hide_colorbar = TRUE,
                 branches_lwd = 0.1,
                 fontsize_row = 5, fontsize_col = 5,
                 heatmap_layers = theme(axis.line=element_blank())
  )
  
  klastry
  p
  return(list(klastry,p))
  
} 