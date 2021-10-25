#' @title crossprodBatch function on GPU
#'
#' @param A rectangular matrices
#' @param D rectangular matrix, columns are diagonals
#' @param invertD set to 1 for C = t(A) D^(-1) A
#' @param Nglobal vector of number of global work items
#' @param Nlocal vector of number of local work items
#' @param NlocalCache elements in local cache
#' 
#' @return  C = A^T A or A^T D A or A^T D^(-1) A, output matrices, stacked row-wise 
#' @useDynLib gpuRandom
#' @export


 
crossprodBatch <- function(C,   # must be batch of square matrices 
                           A,
                           D,
                           invertD,
                           Nglobal, 
                           Nlocal, 
                           NlocalCache,
                           Cstartend, Astartend, Dstartend) {
  
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
  
  
  if((NlocalCache - Nlocal[1]*Nlocal[2])<0){
    warning("a larger NlocalCache required")
  }
  
  crossprodBatchBackend(C,A,D,invertD,Cstartend,Astartend,Dstartend, Nglobal,Nlocal, NlocalCache)
  
  invisible()
  
  
}




