normalization <- function(data, normType) {
  
  if(normType == "RMA"){normData = affy::rma(data)}
  if(normType == "GCRMA"){normData = gcrma(data)}
  if(normType == "MAS"){normData = mas5(data)}

  return(normData)

}

