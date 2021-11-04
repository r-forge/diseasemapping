#' @title createStreamsCpu
#'
#' @description create streams stored on the CPU.
#' 
#' @param n number of streams to create
#' @param initial initial state of first stream, length 6, recycled if shorter
#' @return A stream object on CPU.
#' @useDynLib gpuRandom
#' @export





createStreamsCpu = function(n, 
                            initial){
  
  if(missing(initial)) {
    initial = as.integer(as.integer(2^31-1)*(2*stats::runif(6) - 1) ) 
  }else{
    initial = rep_len(initial, 6)
  }
  
  streamsR <- gpuRandom:::createStreamsCpuBackend(n, initial)
  
  streamsR
  
}

