#' @title Random number generation on GPU
#'
#' @description Generate uniform or Normal random numbers on GPU.
#' 
#' @param x A Vclvector of returned random numbers
#' @param workgroupSize Global dimensions of work-item domain
#' @param localSize Dimensions of work groups
#' @param streams Vectors of integers, 18x workgroupSize
#' @return Altered in-place
#' @useDynLib gpuRandom
#' @export

gpuRn = function(
  x, 
  streams, 
  workgroupSize,   
  localSize, 
  random_type=c("uniform","normal")) {

  verbose = TRUE
  
	if(missing(streams)) {
	  if(missing(workgroupSize)) 
	    {workgroupSize = c(4,64)}
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
  
  if (missing(localSize)){
    localSize = pmax(2,c(localSize, 2, 2)[1:2])}
  
  
  if(verbose) {
    cat('local sizes ', toString(localSize), '\nglobal sizes ', toString(workgroupSize),
        '\n streams ', toString(dim(streams)), '\n')
      }
   
  cpp_gpuRn(x,  streams, workgroupSize, localSize, random_type[1], gpuR::typeof(x))

	invisible(streams)
}














