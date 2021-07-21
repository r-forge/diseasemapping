#' @title qqnorm on GPU
#' 
#' @param y the data sample
#' @param ylim  limits on the plot region
#' @param mu  mean of Normal distribution
#' @param sigma  variance of Normal distribution
#' @param lowertail logical, whether use lower tail probability, default is TRUE
#' @param main plot label
#' @param xlab plot label
#' @param ylab plot label
#' @param Nglobal the size of the index space for use, default is 64 by 4
#' @param Nlocal the work group size of the index space
#' 
#' @return the Normal Q-Q plot
#' 
#' @return the result matrices of LDB in batch 
#'
#' @useDynLib gpuRandom
#' @export

qqnorm<-function(y, ylim, mu=0, sigma=1, lowertail=1,
                  main = "Normal Q-Q Plot",
                  xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
                  Nglobal, Nlocal = c(2, 2),
                  verbose=FALSE, ...){
   
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

    if(missing(Nglobal)) 
    {Nglobal = c(64,4)}
    
    
    if(verbose) {
      cat('local sizes ', toString(Nlocal), '\nglobal sizes ', toString(Nglobal), '\n')
    }
    

 #   p <-gpuR::vclVector(ppoints(n), type=gpuR::typeof(y))
    out <-gpuR::vclVector(length=as.integer(n), type=gpuR::typeof(y))
   
    x <- as.vector(cpp_gpu_qqnorm(out, mu,sigma, lowertail, Nglobal , Nlocal))
    
    x<-x[order(order(as.vector(y)))]  ###
    
    
    if(has.na) {
      y <- x; 
      x <- yN; 
      x[!ina] <- y;
      y <- yN
    }
    
    
      plot(x, y, main = main, xlab = xlab, ylab = ylab, ylim = ylim, ...)
   
       invisible(list(x = x, y = y))
       
      
  }



































