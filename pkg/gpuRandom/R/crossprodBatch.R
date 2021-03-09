#' @title crossprodBatch function on GPU
#'

#' 
#' @return  C = A^T A or A^T D A or A^T D^(-1) A 
#' @useDynLib gpuRandom
#' @export


 
crossprodBatch <- function(C,   # must be a batch of square matrices 
                      A,
                      D,
                      invertD,
                      Nglobal, Nlocal, 
                      NlocalCache,
                      Cstartend,
                      Astartend,
                      Dstartend) {
  
  Nbatches = nrow(C)/ncol(C)
  
  if(missing(Cstartend)) {
    Cstartend=c(0, ncol(C), 0, ncol(C))
  }
  
  if(missing(Astartend)) {
    Astartend=c(0, nrow(A)/Nbatches, 0, ncol(A))
  }
  
  if(missing(Dstartend)) {
    Dstartend=c(0, 1, 0, ncol(D))
  }
  
  
  
  
  
  
  crossprodBatchBackend(C,A,D,invertD,Cstartend,Astartend,Dstartend, Nglobal,Nlocal, NlocalCache)
  
  invisible()
  
  
}




