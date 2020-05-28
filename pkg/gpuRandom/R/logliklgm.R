#' @title Log-likelihood function for Gaussian random fields
#'
#' @useDynLib gpuRandom
#' @export




# 1, \log |\sigma^2 matern() + \tau^2 I|

# 2, (y-D\beta)^\T(\sigma^2V)^(-1)(y-D\beta)


loglikGpu = function(data,  # An object of class SpatialPointsDataFrame
                     paramBatch, # a vclMatrix of parameters shape, range, variance, nugget, anisoRatio, anisoAngleRadians
                     workgroupSize,
                     localSize,
                     verbose=FALSE
                     ){
  
  #coordinates
  
  coordsGpu = vclMatrix(data@coords, nrow(data@coords), ncol(data@coords), type=gpuR::typeof(paramBatch))
  
  #results
  outputBatch = vclMatrix(0, nrow(paramsGpu)*nrow(coordsGpu), nrow(coordsGpu), type=gpuR::typeof(paramsBatch))
  
  gpuRandom:::maternBatchBackend(outputBatch, coordsGpu, paramBatch,  
                                Nglobal = workgroupSize, Nlocal = localSize)
  
  
  
  
  
  
  
  
  
  
  
  
}















