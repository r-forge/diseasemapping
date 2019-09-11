#' @title gpuFisher_test
#'
#' @description Peform Fisher's simulation test on GPU
#'
#' @param sr a Vclvector
#' @param sc a Vclvector
#' @param x a Vclvector
#' @param streams vector of integers, 18x workgroupSize
#' @param workgroupSize global dimensions of work items
#' @param localSize dimensions of work groups
#' @return altered in-place
#' @useDynLib gpuRandom
#' @export

fisher_simGpu=function(
  sr, #marginal row total 
  sc, #marginal column total
  x, # the vector to store test statistics,//extraX,
  streams, 
  workgroupSize,
  localSize){
  
  verbose = TRUE
  
  
  if(missing(streams)) {
    if(missing(workgroupSize)) 
    {workgroupSize = c(64,4)
    streams = cpp_mrg31k3pCreateStreams(prod(workgroupSize))	}	
  } else {
    if(!is.matrix(streams)) {
      warning("streams should be a matrix")
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
  cpp_gpuFisher_test(sr, sc, x, #extraX, 
                     streams, workgroupSize,localSize)
  
  invisible(streams)
  
}



















