distGpu = function(
    x, 
    param = c(range = 1, variance = 1, shape = 1), 
    result = vclMatrix(
        data=0, 
        nrow(x), nrow(x),
        type='double',
        ctx_id = x@.context_index),	
    type = c("variance", "cholesky", "precision", "inverseCholesky")
) {
  
  maxWorkGroupSize <- switch(
      deviceType(result@.platform_index, result@.device_index),
      "gpu" = gpuInfo(result@.platform_index, result@.device_index)$maxWorkGroupSize,
      "cpu" = cpuInfo(result@.platform_index, result@.device_index)$maxWorkGroupSize,
      stop("unrecognized device type")
  )
  
  param = geostatsp::fillParam(param)
  
  file = '/home/patrick/workspace/diseasemapping/pkg/geostatsgpu/inst/CL/dist.cl'
  if(!file.exists(file))
    file <- system.file("CL", "dist.cl", package = "geostatsgpu")
  
  
  if(!file_test("-f", file)){
    stop("kernel file does not exist")
  }
  kernel <- readChar(file, file.info(file)$size)
  
cpp_distGpu(
      x@address, 
      result@address, 
      param, 
      kernel,
      floor(sqrt(maxWorkGroupSize)), 
      x@.context_index - 1)
  
  result  
}


maternGpuOld = function(
    x, 
    param = c(range = 1, variance = 1, shape = 1), 
    result = vclMatrix(
        data=0, 
        nrow(x), nrow(x),
        type='double',
        ctx_id = x@.context_index),	
    type = c("variance", "cholesky", "precision", "inverseCholesky")
) {
  
  maxWorkGroupSize <- switch(
      deviceType(x@.platform_index, x@.device_index),
      "gpu" = gpuInfo(x@.platform_index, x@.device_index)$maxWorkGroupSize,
      "cpu" = cpuInfo(x@.platform_index, x@.device_index)$maxWorkGroupSize,
      stop("unrecognized device type")
  )
  
  param = geostatsp::fillParam(param)
  
 # .Call("cpp_maternGpu", 
 #     x@address, 
 #     result@address, 
 #     param,
 #     sqrt(maxWorkGroupSize), 
 #     x@.context_index - 1)
  
  result  
}
