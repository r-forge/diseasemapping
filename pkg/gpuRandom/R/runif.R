#' @title runif
#' 
#' @param n a number or a vector specifying the size of output matrix
#' @param streams streams matrix of random number seeds
#' @param Nglobal vector of length 2, number of work items
#' @param type random number type, one of "uniform" or "normal"
#' 
#' @return a vclVector or vclMatrix of uniform random numbers
#' 
#' @useDynLib gpuRandom
#' @export



runif = function(
  n, 
  streams, 
  Nglobal,
  type=c("double","float","integer")) {
  
  
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
    if(missing(Nglobal)) {
      Nglobal = c(64,8)
      seedR = as.integer(as.integer(2^31-1)*(2*stats::runif(6) - 1) ) 
      seed <- gpuR::vclVector(seedR, type="integer")  
      streams<-vclMatrix(0L, nrow=512, ncol=12, type="integer")
      CreateStreamsGpuBackend(seed, streams, keepInitial=1)
    }else{
      seedR = as.integer(as.integer(2^31-1)*(2*stats::runif(6) - 1) ) 
      seed <- gpuR::vclVector(seedR, type="integer")  
      streams<-vclMatrix(0L, nrow=prod(Nglobal), ncol=12, type="integer")
      CreateStreamsGpuBackend(seed, streams, keepInitial=1)
    }
  }else if(missing(Nglobal)){
    stop("number of work items needs to be same as number of streams")
  }else if(prod(Nglobal) != nrow(streams)){
      warning("number of work items needs to be same as number of streams")
  }
  


  xVcl<-gpuR::vclMatrix(0L, nrow=n[1], ncol=n[2], type=type[1])    
  
  
  gpuRnBackend(xVcl,streams,Nglobal,"uniform") 
  
  invisible(streams)
  
  if(ncol(xVcl)==1) {
    if (type == "float"){
      invisible(capture.output(xVcl <- as.vclVector(xVcl)))     # an unneeded message from gpuR if float
    }else{
      xVcl <- as.vclVector(xVcl)
    }
      
  }
  
  xVcl
  
}














  