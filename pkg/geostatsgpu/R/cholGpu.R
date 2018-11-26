# To DO:
# find multiple, shared memory in R
# set device sizes from R
# deal with shared memory too large
# split dLocal?
cholGpu = function(x,D=NULL) {

	if(is.null(D)) {
		D = vclVector(0, nrow(x),
		    type='double',
    		ctx_id = x@.context_index)
	}

	maxWorkGroupSize <- switch(
  deviceType(x@.platform_index, x@.device_index),
  "gpu" = gpuInfo(x@.platform_index, x@.device_index)$maxWorkGroupSize,
  "cpu" = cpuInfo(x@.platform_index, x@.device_index)$maxWorkGroupSize,
  stop("unrecognized device type")
  )

maxWorkGroupSize = min(c(maxWorkGroupSize, nrow(x)))

file <- system.file("CL", "cholGpu.cl", package = "geostatsgpu")
kernel <- readChar(file, file.info(file)$size)

	fromC = cpp_cholGpu(
		x@address, 
  		D@address,
  		maxWorkGroupSize, 
		x@.context_index - 1,
		kernel
		)

	res = list(L = x, D=D, logDet = fromC)
	res
}