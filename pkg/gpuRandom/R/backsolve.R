#' @title Backsolve function on GPU
#'
#' @param A a vclMatrix 
#' @param B a vclMatrix 
#' @param diagIsOne logical, indicates if all the diagonal entries in A are one
#' 
#' @return a vclMatrix C that satisfies A*C=B
#' @useDynLib gpuRandom
#' @export



backsolveBatch <- function(C, #output vclmatrix
                           A, B,  #vclmatrices
                           Cstartend,
                           Astartend=c(0, ncol(A), 0, ncol(A)),
                           Bstartend,
                           numbatchB,
                           diagIsOne,
                           workgroupSize, localsize,
                           NlocalCache,
                           verbose=FALSE){


  if(missing(workgroupSize)) {
    workgroupSize <- c(64,64,16)
  }
  
  if(missing(localSize)) {
    localSize <- c(4,4,4)
  }

  if(missing(NlocalCache)) {
    NlocalCache <- 18
  }
  
  
  
  C = vclMatrix(data=0, ncol(A), ncol(B), type=gpuR::typeof(A))
  
  
  
  
  
  
  if(verbose){ message(paste('global work items', workgroupSize,
                             'local work items', localSize))}


 

  backsolveBatchBackend(C, A, B, diagIsOne, workgroupSize, localsize, NlocalCache)


  C

  }

  