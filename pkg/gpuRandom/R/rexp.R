#' @title rexp
#' 
#' @param n a number or a vector specifying the size of output matrix
#' @param streams streams 
#' @param rate mean equals to 1/rate
#' @param Nglobal vector of length 2, number of work items
#' @param type random number type, one of "uniform" or "normal"
#' 
#' @return a vclVector or vclMatrix of exponential random numbers
#' 
#' @useDynLib gpuRandom
#' @export




rexp = function(
  n, 
  rate=1,
  streams, 
  Nglobal,
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

  if(Nglobal[2]<2){
    stop("Nglobal[2] should be larger than 1")
  }

  if(rate <= 0 || !is.finite(rate)){
    stop("invalid rate value")
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
  
  
  
  
  xVcl<-gpuR::vclMatrix(0, nrow=n[1], ncol=n[2], type=type[1])     
  
  
  gpuRnBackend(xVcl,streams, Nglobal,"exponential") 
  
  invisible(streams)
  
  if(ncol(xVcl)==1) xVcl = xVcl[,1]
  
  rate * xVcl
  
}

