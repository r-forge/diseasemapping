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
  
  maxWorkGroupSize <- switch(
      deviceType(x@.platform_index, x@.device_index),
      "gpu" = gpuInfo(x@.platform_index, x@.device_index)$maxWorkGroupSize,
      "cpu" = cpuInfo(x@.platform_index, x@.device_index)$maxWorkGroupSize,
      stop("unrecognized device type")
  )
  
  param = geostatsp::fillParam(param)
  
  .Call("cpp_maternGpu", 
      x@address, 
      result@address, 
      param,
      maxWorkGroupSize, 
      x@.context_index - 1)
  
  result  
}
