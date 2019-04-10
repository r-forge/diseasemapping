#' @title runif gpu
#'
#' @description Random uniforms
#'
#' @param x vector
#' @param workgroupSize workgroup size
#' @param localSize local size
#' @param streams vector of integers, 18x workgroupSize
#' @return altered in-place
#' @useDynLib clRNG
#' @export

rnumberGpu = function(
  x, 
  streams, 
  workgroupSize=1, 
  localSize=1, 
  random_type=c("uniform","normal")) {

	if(missing(streams)) {
		streams = cpp_mrg31k3pCreateStreams(workgroupSize)			
	} else {
		if(!is.matrix(streams)) {
			warning("streams should be a matrix")
		}
		if(ncol(streams) != 18) {
			warning("streams needs 18 columns")
		}
		workgroupSize = ncol(streams)
		# make a deep copy
		streams = matrix(as.vector(streams), nrow(streams), ncol(streams), FALSE, dimnames(streams))
	}
   
	cpp_random_numberGpu(x,  streams, as.integer(workgroupSize), as.integer(localSize), random_type, gpuR::typeof(x))
   ### why?

	invisible(streams)
}














