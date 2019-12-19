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
  
  x<-matrix(0,nrow=n[1],ncol=n[2])
  
  
  if(missing(streams)) {
    if(missing(workgroupSize)) {
      workgroupSize = c(64,8)
      streams = cpp_mrg31k3pCreateStreams(prod(workgroupSize))	
    }else{
      streams = cpp_mrg31k3pCreateStreams(prod(workgroupSize)) 
    }
  }else {
    if(!is.matrix(streams)) {
      warning("streams should be a matrix") }
    
    if(prod(workgroupSize) != nrow(streams))
      warning("number of work items needs to be same as number of streams")
    # make a deep copy
    streams = matrix(as.vector(streams), nrow(streams), ncol(streams), FALSE, dimnames(streams))
  }

    xVcl<-gpuR::vclMatrix(x,type=type[1])###change
  
  
    gpuRnBackend(xVcl,streams,workgroupSize,"normal") 
    
    invisible(streams)
    
    if(ncol(xVcl)==1) xVcl = xVcl[,1]
  xvcl

}
