simEpidemic <- function(params, N, days=10,
  probOnsetMissing=0.7) {
  
  if(!any(names(params[["OnsMedM"]])=="scale"))
    warning("can't find scale parameter, use addScaleParameters")
  
  colNames <- c("infect", "onset",  "med", "hospIn",
     "hospOut", "death", "censor", "type","observedType")
  result <- matrix(NA, N, length(colNames), 
    dimnames = list(NULL, colNames))
  result <- as.data.frame(result)
  
  result$med = sample(days, N, replace=T)
  theTypes = c("M","S","D")
  result$type = sample(factor(theTypes, levels=theTypes),
    N, replace=T, prob=params$probs)
  
  # onset times
  haveOnset = as.logical(rbinom(N, 1, 1-probOnsetMissing))
  
  result[haveOnset,"onset"]  = - round(rweibull(sum(haveOnset),
    shape= getVecParams(params, "OnsMed", "shape")[as.character(result$type)],
    scale= getVecParams(params, "OnsMed", "scale")[as.character(result$type)]
    ) )   
  
    
  result



}


