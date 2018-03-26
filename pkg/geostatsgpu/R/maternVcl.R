#' @title Matern correlation
#'
#' @description Compute Matern covariance function
#'
#' @param x points
#' @param param parameters
#' @param result VCL matrix to store result
#' @return a VCL matrix
#' @examples
#' x = as.matrix(expand.grid(seq(1,by=1,len=23), seq(201,by=1,len=14)))
#' coordsV = gpuR::vclMatrix(x)
#' D3 <- vclMatrix(
#'     data=-1, 
#'     nrow(coordsV), nrow(coordsV),
#'     type='double'
#' )
#' 
#' myParams = c(shape=4.5, range=1.5, variance = 2, nugget = 0, anisoRatio = 2, anisoAngleRadians = pi/4)

#' system.time(v1 <- maternGpu(coordsV, myParams, D3, 'variance'))
#' system.time(mmat <- geostatsp::matern(sp::SpatialPoints(as.matrix(coordsV)), myParams))
#' as.matrix(v1)[1:5,1:5]
#' mmat[1:5,1:5]
#' 
#' @export
maternGpu = function(
    x, 
    param = c(range = 1, variance = 1, shape = 1), 
    result = vclMatrix(
        data=0, 
        nrow(x), nrow(x),
        type='double',
        ctx_id = x@.context_index),	
    type = c("variance", "cholesky", "precision", "inverseCholesky")
) {
  
  type = gsub("iance$|esky$|ision", "", tolower(type)[1])    
  type = c(var=1,chol=2,prec=3,inversechol=4)[type]    

  maxWorkGroupSize <- switch(
      deviceType(result@.platform_index, result@.device_index),
      "gpu" = gpuInfo(result@.platform_index, result@.device_index)$maxWorkGroupSize,
      "cpu" = cpuInfo(result@.platform_index, result@.device_index)$maxWorkGroupSize,
      stop("unrecognized device type")
  )
  maxWorkGroupSize = floor(sqrt(maxWorkGroupSize))
  
  param = geostatsp::fillParam(param)
  
    file <- system.file("CL", "matern.cl", package = "geostatsgpu")

  kernel <- readChar(file, file.info(file)$size)

# use https://github.com/markholland/cholesky/blob/master/OpenCl/kernels.c instead?
  cholFile <- system.file("CL", "dcholesky.cl", package = "gpuR")
  cholKernel <- readChar(cholFile, file.info(cholFile)$size)
  
  upper = FALSE
  cpp_maternGpu(
      x@address, 
      result@address, 
      param,
      type,
      upper,
      kernel,
      cholKernel,
      maxWorkGroupSize, 
      x@.context_index - 1)
  
  if(type == Inf) {
    # TO DO: pass type to cpp_maternGpu
    
    gpuR:::cpp_vclMatrix_custom_chol(
        result@address,
        TRUE, # vcl
        FALSE, # upper
        cholKernel,
        maxWorkGroupSize,
        8L,
        result@.context_index - 1)
  }
  
  result  
}

maternGpuSI = function(
    x, 
    param = c(range = 1, variance = 1, shape = 1), 
    result = vclMatrix(
        data=0, 
        nrow(x), nrow(x),
        type='double',
        ctx_id = x@.context_index), 
    type = c("variance", "cholesky", "precision", "inverseCholesky")
) {
  
  type = gsub("iance$|esky$|ision", "", tolower(type)[1])    
  type = c(var=1,chol=2,prec=3,inversechol=4)[type]
  print(type)
  maxWorkGroupSize <- switch(
      deviceType(result@.platform_index, result@.device_index),
      "gpu" = gpuInfo(result@.platform_index, result@.device_index)$maxWorkGroupSize,
      "cpu" = cpuInfo(result@.platform_index, result@.device_index)$maxWorkGroupSize,
      stop("unrecognized device type")
  )
  
  param = geostatsp::fillParam(param)
  
    file <- system.file("CL", "maternSingleIndex.cl", package = "geostatsgpu")
  kernel <- readChar(file, file.info(file)$size)


  junkD = junkY = junkC = gpuR::vclMatrix(matrix(1.1, 2, 2))

  upper = FALSE
  cpp_maternGpuSingleIndex(
      result@address, 
      junkD@address, junkY@address, junkC@address,
      x@address, 
      param,
      type,
      upper,
      kernel,
      maxWorkGroupSize, 
      x@.context_index - 1)
  
  result  
}
