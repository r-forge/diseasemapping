#' @title runif
#' @useDynLib gpuRandom
#' @export


runif = function(
  n, 
  streams, 
  workgroupSize,
  type=c("double","float")) {
  
  
  if(length(n)>=3){
    stop("'n' has to be a vector of no more than two elements")
  }
  if(length(n)==0){
    stop("specify the number of rows and columns of the output matrix")
  }
  if(length(n)==1){
    n<-c(n,1)
  }
  

  
  
  if(missing(streams)) {
    if(missing(workgroupSize)) {
      workgroupSize = c(64,8)
      result = matrix(0, nrow=prod(workgroupSize), ncol=18)
      colnames(result) = c("current.g1.1", "current.g1.2", "current.g1.3", "current.g2.1", "current.g2.2", "current.g2.3",
                           "initial.g1.1", "initial.g1.2", "initial.g1.3", "initial.g2.1", "initial.g2.2", "initial.g2.3",
                           "substream.g1.1", "substream.g1.2", "substream.g1.3", "substream.g2.1", "substream.g2.2", "substream.g2.3")
      streams = gpuR::vclMatrix(cpp_mrg31k3pCreateStreams(result))
    }else{
      result = matrix(0, nrow=prod(workgroupSize), ncol=18)
      colnames(result) = c("current.g1.1", "current.g1.2", "current.g1.3", "current.g2.1", "current.g2.2", "current.g2.3",
                           "initial.g1.1", "initial.g1.2", "initial.g1.3", "initial.g2.1", "initial.g2.2", "initial.g2.3",
                           "substream.g1.1", "substream.g1.2", "substream.g1.3", "substream.g2.1", "substream.g2.2", "substream.g2.3")
      streams = gpuR::vclMatrix(cpp_mrg31k3pCreateStreams(result))
    }
  }else {
    if(!isS4(streams)) {
      warning("streams should be a S4 matrix") }
    
    if(prod(workgroupSize) != nrow(streams))
      warning("number of work items needs to be same as number of streams")
    # make a deep copy
   # streams = gpuR::vclMatrix(as.matrix(streams), nrow(streams), ncol(streams), FALSE, dimnames(streams))
  }
  
  xVcl<-gpuR::vclMatrix(0, nrow=n[1], ncol=n[2], type=type[1])    
  
  
  gpuRnBackend(xVcl,streams,workgroupSize,"uniform") 
  
  invisible(streams)
  
  if(ncol(xVcl)==1) xVcl <- as.vclVector(xVcl)
  xVcl
  
}
