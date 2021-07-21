#' @title multiplyLowerDiagonalBatch function on GPU
#'
#' @param output the output = LDB
#' @param L  lower triangular matrices in batce
#' @param D  diagonal matrices in batch
#' @param B  matrices in batch
#' @param diagIsOne logical, whether the diagonal of L is one 
#' @param transformD how to transform D
#' @param Nglobal the size of the index space for use
#' @param Nlocal the work group size of the index space 
#' @param NlocalCache a number
#' 
#' @return the result matrices of LDB in batch 
#' 
#' @useDynLib gpuRandom
#' @export




# output = L  D B, L lower triangular, D diagonal

multiplyLowerDiagonalBatch <- function(
                      output, L, D, B,
                      diagIsOne, # diagonal of L is one
                      transformD, 
                      Nglobal,
                      Nlocal,
                      NlocalCache){
  
  
  
  gpuRandom:::multiplyLowerDiagonalBatchBackend(
               output,
               L,
               D,
               B,
               diagIsOne,    
               transformD,
               Nglobal,
               Nlocal,
               NlocalCache)
  
  
  
}











