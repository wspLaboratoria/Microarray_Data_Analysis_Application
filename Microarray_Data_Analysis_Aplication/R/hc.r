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
  
  mydatascale <- (scale(t(dataHeatmap))) # Centers
  d <- dist(mydatascale, method = dist_method)
  hc <- hclust(d, method =clust_method)
  
  sub_grp <- cutree(hc, k = cut_number)
  
  klastry <- fviz_cluster(list(data = mydatascale, cluster = sub_grp), geom='text')

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
  
  p
  
  return(list(klastry,p))

}