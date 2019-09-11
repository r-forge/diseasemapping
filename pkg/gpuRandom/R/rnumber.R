#' @title runif gpu
#'
#' @description Random uniforms
#'
#' @param x vector
#' @param workgroupSize workgroup size
#' @param localSize local size
#' @param streams vector of integers, 18x workgroupSize
#' @return altered in-place
#' @export

rnumberGpu = function(
  x, 
  streams, 
  workgroupSize,   
  localSize=c(2,2), 
  random_type=c("uniform","normal")) {

  verbose = TRUE
  
	if(missing(streams)) {
	  if(missing(workgroupSize)) 
	    {workgroupSize = c(2,64)}
	  workgroupSize = pmax(2, c(workgroupSize, 2, 2)[1:2])
	  streams = cpp_mrg31k3pCreateStreams(prod(workgroupSize))			
	} else {
		if(!is.matrix(streams)) {
			warning("streams should be a matrix")
		}
		if(ncol(streams) != 18) {
			warning("streams needs 18 columns")
		}
	  if(missing(workgroupSize)) {
	    workgroupSize = c(nrow(streams)/2,2)
	  }
		if(prod(workgroupSize) != nrow(streams))
		  warning("number of work items needs to be same as number of streams")
		# make a deep copy
		streams = matrix(as.vector(streams), nrow(streams), ncol(streams), FALSE, dimnames(streams))

			}
  
  
  localSize = pmax(2,c(localSize, 2, 2)[1:2])
  
  if(verbose) {
    cat('local sizes ', toString(localSize), '\nglobal sizes ', toString(workgroupSize),
        '\n streams ', toString(dim(streams)), '\n')
      }
   
	cpp_random_numberGpu(x,  
	                     streams, 
	                     workgroupSize,
                       localSize,
	                     random_type[1], 
	                     gpuR::typeof(x))

	invisible(streams)
}














