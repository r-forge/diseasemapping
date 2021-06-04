#' @title Create matern Params on GPU
#'
#' @useDynLib gpuRandom
#' @export



maternParamsGpu = function(
  myParamsBatch, type = c('float','double')){
  

  myParamsBatch = t(apply(myParamsBatch, 1, geostatsp::fillParam))
  myParamsBatch = cbind(myParamsBatch, matrix(0, nrow(myParamsBatch), 22-ncol(myParamsBatch)))
  paramsGpu =  vclMatrix(myParamsBatch, type=theType)
  
  paramsGpu
  
}







