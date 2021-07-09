#' @title multiplyLowerDiagonalBatch function on GPU
#'
#' @useDynLib gpuRandom
#' @export




# output = L  D B, L lower triangular, D diagonal

multiplyLowerDiagonalBatch <- function(
                      output, L, D, B,
                      diagIsOne, # diagonal of L is one
                      transformD, 
                      workgroupSize,
                      localSize,
                      NlocalCache){
  
  
  
  gpuRandom:::multiplyLowerDiagonalBatchBackend(
               output,
               L,
               D,
               B,
               diagIsOne,    
               transformD,
               workgroupSize,
               localSize,
               NlocalCache)
  
  
  
}











