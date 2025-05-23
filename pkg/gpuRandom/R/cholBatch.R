#' @title Backsolve function on GPU
#'
#' @param A a vclMatrix, positive definite
#' @param D a vclMatrix, diagonal elements of LDL^t 
#' @param numbatchD number of batches in D
#' 
#' @return  vclMatrix L and D that satisfies A=LDL^T
#' @useDynLib gpuRandom
#' @export



cholBatch <- function(A,
                      D,
                      numbatchD,
                      Nglobal,
                      Nlocal,   # needs Nglobal[2]=Nlocal[2]
                      NlocalCache,
                      Astartend,
                      Dstartend,
                      verbose=FALSE){
  

  
  if(missing(Astartend)) {
    Astartend=c(0, nrow(A)/numbatchD, 0, ncol(A))
  }
  
  if(missing(Dstartend)) {
    Dstartend=c(0, numbatchD, 0, ncol(D))
  }
  
  
   if(Nlocal[2]!=Nglobal[2]){
     warning("local and global work sizes should be identical for dimension 2, ignoring global")
     Nglobal[2]=Nlocal[2]
   }
  
 
  
  if(verbose){ message(paste('global work items', Nglobal,
                             'local work items', Nlocal))}
  
  
  
  
  cholBatchBackend(A, D, 
                   Astartend, Dstartend, 
                   numbatchD,
                   Nglobal, Nlocal, NlocalCache)
  
   #theResult = list(L=A, diag=D)


   #theResult
   invisible()
  
  
}