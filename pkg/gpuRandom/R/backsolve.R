#' @title Backsolve function on GPU
#'
#' @param C a vclMatrix
#' @param A a vclMatrix, unit lower triangular
#' @param B a vclMatrix 
#' @param numbatchB number of batches in B
#' @param diagIsOne logical, indicates if all the diagonal entries in A are one
#' @return a vclMatrix C that satisfies A*C=B
#' @useDynLib gpuRandom
#' @export



backsolveBatch <- function(C, 
                           A,  # must be batches of square matrices
                           B,  #vclmatrices
                           numbatchB, #sometimes B can have only 1 batch, for repeated same batches
                           diagIsOne,
                           Nglobal, 
                           Nlocal, 
                           NlocalCache,
                           Cstartend,
                           Astartend,
                           Bstartend,
                           verbose=FALSE){

  
  nbatch<-nrow(A)/ncol(A)

  if(missing(Cstartend)) {
    Cstartend=c(0, nrow(C)/nbatch, 0, ncol(C))
  }
  
  if(missing(Astartend)) {
    Astartend=c(0, nrow(A)/nbatch, 0, ncol(A))
  }
  
  if(missing(Bstartend)) {
    Bstartend=c(0, nrow(B)/numbatchB, 0, ncol(B))
   }
  

  
  
  
  
  
  if(verbose){ message(paste('global work items', Nglobal,
                             'local work items', Nlocal))}


 

  backsolveBatchBackend(C, A, B, 
                        Cstartend, Astartend, Bstartend, 
                        numbatchB, diagIsOne, 
                        Nglobal, Nlocal, NlocalCache)



  }

  