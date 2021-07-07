#' @title rnorm
#' @useDynLib gpuRandom
#' @export




rnorm = function(
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
  
  #x<-matrix(0,nrow=n[1],ncol=n[2])
  
  if(missing(streams)) {
    if(missing(workgroupSize)) {
      workgroupSize = c(64,8)
      seedR = as.integer(as.integer(2^31-1)*(2*stats::runif(6) - 1) ) 
      seed <- gpuR::vclVector(seedR, type="integer")  
      streams<-vclMatrix(0L, nrow=512, ncol=18, type="integer")
      CreateStreamsGpuBackend(seed, streams, keepInitial=1)
    }else{
      seedR = as.integer(as.integer(2^31-1)*(2*stats::runif(6) - 1) ) 
      seed <- gpuR::vclVector(seedR, type="integer")  
      streams<-vclMatrix(0L, nrow=prod(workgroupSize), ncol=18, type="integer")
      CreateStreamsGpuBackend(seed, streams, keepInitial=1)
    }
  }else if(missing(workgroupSize)){
    stop("number of work items needs to be same as number of streams")
  }else if(prod(workgroupSize) != nrow(streams)){
    warning("number of work items needs to be same as number of streams")
  }
  
  

  
  xVcl<-gpuR::vclMatrix(0, nrow=n[1], ncol=n[2], type=type[1])     
  
  
  gpuRnBackend(xVcl,streams,workgroupSize,"normal") 
  
  invisible(streams)
  
  if(ncol(xVcl)==1) xVcl = xVcl[,1]
  
  xVcl
  
}








# if(missing(streams)) {
#   if(missing(workgroupSize)) {
#     workgroupSize = c(64,8)
#     result = matrix(0L, nrow=prod(workgroupSize), ncol=18)
#     colnames(result) = c("current.g1.1", "current.g1.2", "current.g1.3", "current.g2.1", "current.g2.2", "current.g2.3",
#                          "initial.g1.1", "initial.g1.2", "initial.g1.3", "initial.g2.1", "initial.g2.2", "initial.g2.3",
#                          "substream.g1.1", "substream.g1.2", "substream.g1.3", "substream.g2.1", "substream.g2.2", "substream.g2.3")
#     streams = gpuR::vclMatrix(cpp_mrg31k3pCreateStreams(result))
#   }else{
#     result = matrix(0L, nrow=prod(workgroupSize), ncol=18)
#     colnames(result) = c("current.g1.1", "current.g1.2", "current.g1.3", "current.g2.1", "current.g2.2", "current.g2.3",
#                          "initial.g1.1", "initial.g1.2", "initial.g1.3", "initial.g2.1", "initial.g2.2", "initial.g2.3",
#                          "substream.g1.1", "substream.g1.2", "substream.g1.3", "substream.g2.1", "substream.g2.2", "substream.g2.3")
#     streams = gpuR::vclMatrix(cpp_mrg31k3pCreateStreams(result))
#   }
# }else if(missing(workgroupSize)){
#   stop("number of work items needs to be same as number of streams")
# }else{
#   if(!isS4(streams)) {
#     warning("streams should be a S4 matrix") }
#   
#   if(prod(workgroupSize) != nrow(streams))
#     warning("number of work items needs to be same as number of streams")
#   # make a deep copy
#   # streams = gpuR::vclMatrix(as.matrix(streams), nrow(streams), ncol(streams), FALSE, dimnames(streams))
# }
