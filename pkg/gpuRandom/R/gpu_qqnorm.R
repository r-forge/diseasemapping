#' @title gpu_qqnorm
#'
#' @useDynLib gpuRandom
#' @export

qqnorm<-function(y, ylim, mu, sigma, lowertail=1,
                  main = "Normal Q-Q Plot",
                  xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
                  workgroupSize, localSize=c(2,2),
                  verbose=FALSE){
   
   if(has.na <- any(ina <- is.na(y))) { ## keep NA's in proper places
      yN <- y
      y <- y[!ina]
    }
    if(0 == (n <- length(y)))
      stop("y is empty or has only NAs")
    
    if (missing(ylim))
      ylim <- range(y)
    
    if(sigma<0){
      stop("sigma must not be less than 0")
    }
    
    if(sigma==0){
      x<- rep(mu, n)
    }

    if(missing(workgroupSize)) 
    {workgroupSize = c(64,4)}
    
    
    if(verbose) {
      cat('local sizes ', toString(localSize), '\nglobal sizes ', toString(workgroupSize), '\n')
    }
    

 #   p <-gpuR::vclVector(ppoints(n), type=gpuR::typeof(y))
    out <-gpuR::vclVector(length=as.integer(n), type=gpuR::typeof(y))
   
    x <- as.vector(cpp_gpu_qqnorm(out, mu,sigma, lowertail, workgroupSize , localSize))
    
    x<-x[order(order(as.vector(y)))]
    
    
    if(has.na) {
      y <- x; 
      x <- yN; 
      x[!ina] <- y;
      y <- yN
    }
    
    
      plot(x, y, main = main, xlab = xlab, ylab = ylab, ylim = ylim, ...)
   
       invisible(list(x = x, y = y))
       
      
  }



































