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
  
  print(type)
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
  
  file = '/home/patrick/workspace/diseasemapping/pkg/geostatsgpu/inst/CL/matern.cl'
  if(!file.exists(file))
    file <- system.file("CL", "matern.cl", package = "geostatsgpu")
  
  
  if(!file_test("-f", file)){
    stop("kernel file does not exist")
  }
  kernel <- readChar(file, file.info(file)$size)

#  cholFile <- system.file("CL", "dcholesky.cl", package = "gpuR")
#  cholKernel <- readChar(cholFile, file.info(cholFile)$size)
  
  cpp_maternGpu(
      x@address, 
      result@address, 
      param, 
      kernel,
#      cholKernel,
#      type, 
      floor(sqrt(maxWorkGroupSize)), 
      x@.context_index - 1)
  if(type == 2) {
    .Call('_gpuR_cpp_vclMatrix_custom_chol', PACKAGE = 'gpuR', 
        result, 
        TRUE, TRUE, cholKernel, 
        floor(sqrt(maxWorkGroupSize)), type, x@.context_index - 1)
  }
  
  result  
}

