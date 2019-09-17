#' @title Fisher's exact test on GPU
#'
#' @description Peform Fisher's exact test with Monte Carlo simulation on GPU.
#'
#' @param sr A Vclvector specifying the marginal row totals
#' @param sc A Vclvector specifying the marginal column totals
#' @param x A Vclvector for putting test statistics.
#' @param streams Vectors of integers, 18x workgroupSize
#' @param workgroupSize Global dimensions of work-item domain
#' @return Altered in-place
#' @useDynLib gpuRandom
#' @export

gpuFisher_test=function(
  sr, #marginal row total 
  sc, #marginal column total
  x, #extraX,
  streams, 
  workgroupSize,
  localSize=c(1,1)){
  
  verbose = TRUE
  
  if(missing(streams)) {
    if(missing(workgroupSize)) 
    {workgroupSize = c(64,4)
    streams = cpp_mrg31k3pCreateStreams(prod(workgroupSize))}
  }else{
    if(!is.matrix(streams))
      {warning("streams should be a matrix")}
    
    if(prod(workgroupSize) != nrow(streams))
      warning("number of work items needs to be same as number of streams")

    streams = matrix(as.vector(streams), nrow(streams), ncol(streams), FALSE, dimnames(streams))
     }
  
 
  localSize = pmax(2,c(localSize, 2, 2)[1:2])
  
  if(verbose) {
    cat('local sizes ', toString(localSize), '\nglobal sizes ', toString(workgroupSize),
        '\n streams ', toString(dim(streams)), '\n')
  }
  cpp_gpuFisher_test(sr, sc, x, streams, workgroupSize, localSize)
  
  invisible(streams)
  
}



















