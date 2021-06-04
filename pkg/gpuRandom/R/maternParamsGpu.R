#' @title Create matern Params on GPU
#'
#' @useDynLib gpuRandom
#' @export



maternGpuParam = function(x, type='double') {
  x = geostatsp::fillParam(x)
  paramsGpu = vclMatrix(cbind(x, 
                                matrix(0, nrow(x), 22-ncol(x))),type=type)
  paramsGpu

  
}







