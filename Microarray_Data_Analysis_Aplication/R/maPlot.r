maPlot <- function(data,sample_name){
  
  data.matrix=exprs(data)
  Names=colnames(data.matrix)
  indx=which(Names==sample_name)
  
  
  if (class(data) == 'ExpressionSet') {
    plot <- oligo:: MAplot(data, which=indx)
  }
  if (class(data) == 'AffyBatch') {
    plot <- affy::MAplot(data, which=sample_name)
  }
  return(plot)
}
