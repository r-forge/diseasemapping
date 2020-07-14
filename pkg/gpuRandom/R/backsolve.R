#' @title Backsolve function on GPU
#'
#' @param A a vclMatrix 
#' @param B a vclMatrix 
#' @param diagIsOne logical, indicates if all the diagonal entries in A are one
#' 
#' @return a vclMatrix C that satisfies A*C=B
#' @useDynLib gpuRandom
#' @export



backsolveBatch <- function(A, B,  #vclmatrices
                           diagIsOne, #
                           workgroupSize, localsize,
                           NlocalCache,
                           verbose=FALSE){
  

  if(missing(workgroupSize)) {
    workgroupSize <- c(64,64)
    )



  if(verbose){ message(paste('global work items', workgroupSize, 
                           'local work items', localSize))}


  C = vclMatrix(data=0, ncol(A), ncol(B), type=gpuR::typeof(A))

  backsolveBatchBackend(C, A, B, diagIsOne, workgroupSize, localsize, NlocalCache) 


  C

  }
  
  