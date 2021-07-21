#' @title logfactsumBackend on GPU
#' 
#' 
#' @param x an R matrix
#' @param Nglobal the size of the index space for use
#' 
#' @return   the sum of log-factorials of the matrix
#'
#' @useDynLib gpuRandom
#' @export




logfactSum <- function(x,      # an R matrix
                       Nglobal) {
  
  
  if(any(dim(x) < 2L))
   stop("table must have at least 2 rows and columns")
  
  if(!is.integer(x)) {
          xo <- x
          x <- round(x)
              if(any(x > .Machine$integer.max))
              stop("'x' has entries too large to be integer")
          if(!identical(TRUE, (ax <- all.equal(xo, x))))
           warning(gettext("matrix has been rounded to integer", ax), domain = NA)
           storage.mode(x) <- "integer"
     }
    
    
  x <- gpuR::vclMatrix(x,type="integer")
  
  
  result <- logfactsumBackend(x, Nglobal)
  
  
  result
  
  
  
}




