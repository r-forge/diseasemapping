#' @title Fisher's exact test on GPU
#'
#' @useDynLib gpuRandom
#' @export



fisher.sim=function(
  x,
  B,# number of simualtion,
  type = c('float','double')[1+gpuR::gpuInfo()$double_support],
  returnStatistics = FALSE,
  streams, 
  workgroupSize,
  localSize=c(2,2),
  verbose = FALSE){

  METHOD <- "Fisher's Exact Test for Count Data"
  
  if (is.data.frame(x))
    x <- as.matrix(x)
  
  if(is.matrix(x)) {
    if(any(dim(x) < 2L))
      stop("'x' must have at least 2 rows and columns")
    
    if(!is.numeric(x) || any(x < 0) || any(is.na(x)))
      stop("all entries of 'x' must be nonnegative and finite")
    
    if(!is.integer(x)) {
        xo <- x
        x <- round(x)
        if(any(x > .Machine$integer.max))
        stop("'x' has entries too large to be integer")
        if(!identical(TRUE, (ax <- all.equal(xo, x))))
        warning(gettextf("'x' has been rounded to integer: %d", ax), domain = NA)
        storage.mode(x) <- "integer"
    }
  
   }else{
    stop("'x'  must be matrix")
   }
  
  STATISTIC <- -sum(lfactorial(x))
  ##STATISTIC is negative
  almost.1 <- 1 + 64 * .Machine$double.eps
  
  threshold = STATISTIC/almost.1
  
  if(missing(streams)) {
       if(missing(workgroupSize)) 
           {workgroupSize = c(64,4)
            streams = cpp_mrg31k3pCreateStreams(prod(workgroupSize))	}	
        }else {
              if(!is.matrix(streams)) {
               warning("streams should be a matrix") }
    
    if(prod(workgroupSize) != nrow(streams))
      warning("number of work items needs to be same as number of streams")
    # make a deep copy
    streams = matrix(as.vector(streams), nrow(streams), ncol(streams), FALSE, dimnames(streams))
  }
  
  localSize = pmax(2,c(localSize, 2, 2)[1:2])
  
  if(verbose) {
    cat('local sizes ', toString(localSize), '\nglobal sizes ', toString(workgroupSize),
        '\n streams ', toString(dim(streams)), '\n')
  }
  

  sr0 <- rowSums(x)
  sc0 <- colSums(x)
  ## we drop all-zero rows and columns
  x <- x[sr0 > 0, sc0 > 0, drop = FALSE]

  x<-gpuR::vclMatrix(x, type='int')

  
  print(class(x))
  
  results <- gpuR::vclVector(length=B, type=type)

  PVAL <- NULL
  
  po<-cpp_gpuFisher_test(x, results, threshold, streams, workgroupSize,localSize)
  
  PVAL <- (1 + po ) / (B + 1)
  

  
  if (returnStatistics){
  
  list(pval = PVAL, sim = results, streams=streams)
    
  }else {
    list(pval = PVAL)
    }

  invisible(streams)
  
  
}




















