#' @title Gemmbatch on GPU
#'
#' @param A a vclMatrix consists of batches of submatrices
#' @param B a vclMatrix consists of batches of submatrices
#' @param rowbatch number of batches in row
#' @param Acolbatch number of batches in the column of A
#' @param Bcolbatch number of batches in the column of B
#' 
#' @return a vclMatrix of submatrices multiplication results
#' @useDynLib gpuRandom
#' @export



gemmBatch <- function(
  A, B,  #vclmatrices
  Arowbatch, Browbatch,
  Acolbatch, Bcolbatch,
  need_transpose,
  workgroupSize,
  verbose=FALSE){
  
  
  if( (Arowbatch != Browbatch) & (Browbatch != 1)) 
    stop("A and B must have same number of blocks in row or B can have only one block in row sometimes") 
  
  if( (Arowbatch != Browbatch) & ((Acolbatch != Bcolbatch))) 
    stop("A and B must have same number of row batches or A and B must have same number of col batches") 
  
  
  if((Acolbatch != Bcolbatch) & (Bcolbatch != 1) & (Acolbatch != 1))
    stop("A and B must have same number of blocks in column or either A or B can have only one block in column")
  
  if(missing(workgroupSize)) {
    workgroupSize <- c(64,8,8)
  }
  
  localSize = c(1, 1, 1)
  
  if(verbose){ message(paste('global work items', workgroupSize, 
                             'local work items', localSize))}
  
  
  x <- max(c(Acolbatch,Bcolbatch))
  
  if (need_transpose){
  C = vclMatrix(data=0, ncol(A)/Acolbatch*Arowbatch, ncol(B)/Bcolbatch*x, type=gpuR::typeof(A))
  }else{
  C = vclMatrix(data=0, nrow(A), ncol(B)/Bcolbatch*x, type=gpuR::typeof(A))  
  }

  
  gemmBatchBackend(A, B, C,  
                   Arowbatch, Browbatch, Acolbatch, Bcolbatch,
                   need_transpose,  workgroupSize)
  
  
  C
  
  }








