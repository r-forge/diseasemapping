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

runifGpu = function(x, workgroupSize=1, streams, localSize=1) {

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

	cpp_runifGpu(x,  streams, as.integer(workgroupSize), as.integer(localSize), gpuR::typeof(x))

	invisible(streams)
}














