#' @title logfactsumBackend on GPU
#'
#' @useDynLib gpuRandom
#' @export




logfactSum <- function(x, 
                       workgroupSize) {
  
  
  if(any(dim(x) < 2L))
   stop("table must have at least 2 rows and columns")
  
  if(!is.integer(x)) {
          xo <- x
          x <- round(x)
              if(any(x > .Machine$integer.max))
              stop("'x' has entries too large to be integer")
          if(!identical(TRUE, (ax <- all.equal(xo, x))))
           warning(gettextf("'x' has been rounded to integer: %d", ax), domain = NA)
           storage.mode(x) <- "integer"
     }
    
    
  x <- gpuR::vclMatrix(x,type="integer")
  
  
  result <- logfactsumBackend(x, workgroupSize)
  
  
  result
  
  
  
}




