#' @title maternBatch function on GPU
#'
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




