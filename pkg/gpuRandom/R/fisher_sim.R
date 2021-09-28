#' @title Fisher's exact test on GPU
#'
#'
#' @param x a vclMatrix for Fisher's test
#' @param B requested number of simulations.
#' @param streams streams objects. 
#' @param type "double" or "float" of test statistics.
#' @param Nglobal the size of the index space for use.
#' 
#' 
#' @return a list of results
#' 
#' 
#' @useDynLib gpuRandom
#' @export



fisher.sim=function(
  x, # a vclMatrix
  B, # requested number of simualtion,
  streams, 
  type = c('float','double')[1+gpuR::gpuInfo()$double_support],
  returnStatistics = FALSE,
  Nglobal,
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
  # lfactorial(x)=lgamma(x+1)
  # STATISTIC <- -sum(lfactorial(x))          ##STATISTIC is negative
  # almost.1 <- 1 + 64 * .Machine$double.eps
  # 
  # threshold = STATISTIC/almost.1
  
  if(missing(streams)) {
    if(missing(Nglobal)) {
      Nglobal = c(64,16)
      seedR = as.integer(as.integer(2^31-1)*(2*stats::runif(6) - 1) ) 
      seed <- gpuR::vclVector(seedR, type="integer")  
      streams<-vclMatrix(0L, nrow=1024, ncol=12, type="integer")
      CreateStreamsGpuBackend(seed, streams, keepInitial=1)
      
    }else{
      seedR = as.integer(as.integer(2^31-1)*(2*stats::runif(6) - 1) ) 
      seed <- gpuR::vclVector(seedR, type="integer")  
      streams<-vclMatrix(0L, nrow=prod(Nglobal), ncol=12, type="integer")
      CreateStreamsGpuBackend(seed, streams, keepInitial=1)
    }
  }else {
    if(!isS4(streams)) {
      warning("streams should be a S4 matrix") }
    
    if(prod(Nglobal) != nrow(streams))
      warning("number of work items needs to be same as number of streams")
    # make a deep copy
    # streams = gpuR::vclMatrix(as.matrix(streams), nrow(streams), ncol(streams), FALSE, dimnames(streams))
  }
  
  Nlocal = c(1, 1)
  
  if(verbose) {
    cat('local sizes ', toString(Nlocal), '\nglobal sizes ', toString(Nglobal),
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
  
  TotalWorkItems <- prod(Nglobal)
  
  
  simPerItem1<-B %/% TotalWorkItems  
  if(simPerItem1==0){
    warning("Number of simulations should be larger than number of workitems")
  }
  TotalSim1<-simPerItem1 * TotalWorkItems
  remainder1<-B %% TotalWorkItems
  
  
  if(remainder1 >0){
    dim2 <- ceiling(remainder1/128)
    Nglobal2 <- c(128,dim2)
    TotalSim2 <- prod(Nglobal2)
  }else{
    TotalSim2 <- 0  
  }
  
  
  
  if(returnStatistics) {
    results <- gpuR::vclVector(length=as.integer(TotalSim1+TotalSim2), type=type)
  } else {
    results <- gpuR::vclVector(length=as.integer(1), type=type)
  }
  
  
  PVAL <- NULL
  counts1<-cpp_gpuFisher_test(x, results, 0, simPerItem1, streams, Nglobal, Nlocal)
  
  
  
  if(remainder1 >0){
    seedR <- as.integer(as.integer(2^31-1)*(2*stats::runif(6) - 1) ) 
    seed <- gpuR::vclVector(seedR, type="integer")  
    streams2<-vclMatrix(0L, nrow=prod(Nglobal2), ncol=12, type="integer")
    CreateStreamsGpuBackend(seed, streams2, keepInitial=1)  
    counts2<-cpp_gpuFisher_test(x, results, TotalSim1, 1, streams2, Nglobal2, Nlocal)
  }else{
    counts2 <- 0
  }
  
  PVAL <- (1 + counts1 + counts2) / (TotalSim1 + TotalSim2 + 1)
  #counts<-10
  #PVAL<-0.1
  # format(PVAL, digits=5)
  #if(class(PVAL) == 'try-error') {
  #  PVAL = counts
  #}
  
  
  
  if (returnStatistics){   
    theResult = list(p.value = PVAL, simNum=(TotalSim1 + TotalSim2), sim = results, counts=(counts1 + counts2), streams=streams)
    
  }else {
    
    theResult = list(p.value=PVAL, simNum=(TotalSim1 + TotalSim2), counts=(counts1 + counts2), streams=streams)
  }
  
  theResult
  
  
}

















