#' @title Fisher's exact test on GPU
#'
#' @useDynLib gpuRandom
#' @export



fisher.sim=function(
  x, # a vclMatrix
  B, # number of simualtion,
  streams, 
  type = c('float','double')[1+gpuR::gpuInfo()$double_support],
  returnStatistics = FALSE,
  workgroupSize,
  verbose = FALSE){

  # METHOD <- "Fisher's Exact Test for Count Data"
  
  # if (is.data.frame(x))
  #   x <- as.matrix(x)
  # 
  # if(is.matrix(x)) {
  #   if(any(dim(x) < 2L))
  #     stop("'x' must have at least 2 rows and columns")
  #   
  #   if(!is.numeric(x) || any(x < 0) || any(is.na(x)))
  #     stop("all entries of 'x' must be nonnegative and finite")
  #   
  #   if(!is.integer(x)) {
  #       xo <- x
  #       x <- round(x)
  #       if(any(x > .Machine$integer.max))
  #       stop("'x' has entries too large to be integer")
  #       if(!identical(TRUE, (ax <- all.equal(xo, x))))
  #       warning(gettextf("'x' has been rounded to integer: %d", ax), domain = NA)
  #       storage.mode(x) <- "integer"
  #   }
  # 
  #  }else{stop("'x'  must be matrix")}
  

  
  # STATISTIC <- -sum(lfactorial(x))          ##STATISTIC is negative
  # almost.1 <- 1 + 64 * .Machine$double.eps
  # 
  # threshold = STATISTIC/almost.1
  
  if(missing(streams)) {
    if(missing(workgroupSize)) {
      workgroupSize = c(64,16)
      streams = gpuR::vclMatrix(cpp_mrg31k3pCreateStreams(prod(workgroupSize)))
    }else{
      streams = gpuR::vclMatrix(cpp_mrg31k3pCreateStreams(prod(workgroupSize)))
    }
  }else {
    if(!isS4(streams)) {
      warning("streams should be a S4 matrix") }
    
    if(prod(workgroupSize) != nrow(streams))
      warning("number of work items needs to be same as number of streams")
    # make a deep copy
    # streams = gpuR::vclMatrix(as.matrix(streams), nrow(streams), ncol(streams), FALSE, dimnames(streams))
  }
  
  localSize = c(1, 1)
  
  if(verbose) {
    cat('local sizes ', toString(localSize), '\nglobal sizes ', toString(workgroupSize),
        '\n streams ', toString(dim(streams)), '\n')
  }
  

  # sr0 <- rowSums(x)
  # sc0 <- colSums(x)
  # 
  #   
  # ## we drop all-zero rows and columns
  # x <- x[sr0 > 0, sc0 > 0, drop = FALSE]
  # 
  # xVcl<-gpuR::vclMatrix(x, type='integer') 
  
#  print(class(xVcl))

  
  
  if(returnStatistics) {
    results <- gpuR::vclVector(length=as.integer(B), type=type)
  } else {
    results <- gpuR::vclVector(length=as.integer(1), type=type)
  }
  
  
  PVAL <- NULL
  
  counts<-cpp_gpuFisher_test(x, results, as.integer(B), streams, workgroupSize,localSize)
  
  
  # if(verbose)
  #   print(theTime)
  #time 
  
  PVAL <- ((1 + counts ) / (as.integer(B) + 1))
  
  # format(PVAL, digits=5)
  #if(class(PVAL) == 'try-error') {
  #  PVAL = counts
  #}
  

  
  if (returnStatistics){
  
  theResult = list(p.value = PVAL, sim = results, counts=counts, streams=streams)
    
  }else {
    
  theResult = list(p.value=PVAL, counts=counts, streams=streams)
  }

  theResult
  
  
}




















