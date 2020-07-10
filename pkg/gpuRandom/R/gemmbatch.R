#' @title Gemmbatch on GPU
#'
#' @param A a vclMatrix consists of batches of submatrices
#' @param B a vclMatrix consists of batches of submatrices
#' @param rowbatch number of batches in row
#' @param Acolbatch number of batches in the column of A
#' @param Bcolbatch number of batches in the column of B
#' 
#' @return vclMatrix of submatrices multiplication results
#' @useDynLib gpuRandom
#' @export



gemmBatch <- function(
  A, B,  #vclmatrices
  rowbatch, Acolbatch, Bcolbatch,
  need_transpose,
  workgroupSize,
  verbose=FALSE){
  
  if((Acolbatch != Bcolbatch) & (Bcolbatch != 1))
    stop("A and B should have same number of blocks in column or B should have only one block in column")
  
  localSize = c(1, 1, 1)
  
  if(verbose){ message(paste('global work items', workgroupSize, 
                             'local work items', localSize))}
  
  
  
  if (need_transpose){
  C = vclMatrix(data=0, ncol(A)/Acolbatch*rowbatch, ncol(B)/Bcolbatch*Acolbatch, type=gpuR::typeof(A))
  }else{
  C = vclMatrix(data=0, nrow(A), ncol(B)/Bcolbatch*Acolbatch, type=gpuR::typeof(A))  
  }

  
  gemmBatchBackend(A, B, C,  
                   rowbatch, Acolbatch, Bcolbatch,
                   need_transpose,  workgroupSize)
  
  
  C
  
  }








