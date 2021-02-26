#' @title maternBatch function on GPU
#'
#' @param A a vclMatrix, positive definite
#' @param D a vclMatrix, diagonal elements of LDL^t 
#' @param numbatchD number of batches in D
#' 
#' @return  vclMatrix L and D that satisfies A=LDL^T
#' @useDynLib gpuRandom
#' @export


maternBatch <- function(var,  # the output matern matrices
                        coords,
                        param, #22 columns 
                        Nglobal,
                        Nlocal,
                        startrow,   # new added
                        numberofrows){
  
  
  
  if(missing(startrow) | missing(numberofrows)) {
    startrow=0
    numberofrows=nrow(param)
  }

  
  
  maternBatchBackend(var, coords, param,
                     Nglobal, Nlocal,
                     startrow, numberofrows)
  

  invisible()
  
  
}




