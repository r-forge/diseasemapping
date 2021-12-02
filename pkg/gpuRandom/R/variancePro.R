#' @title Estimate Log-likelihood for Gaussian random fields when stderror are given
#' @import data.table
#' @useDynLib gpuRandom
#' @export


######## this function will be incorporated in the other function later i think

likLgm_variancePro <- function(stderror, #a vector  given by the user 
                               Nobs,  # number of observations.
                               Nparam,
                               Ndata,
                               BoxCox,
                               detVar, # vclVector
                               ssqResidual,   # vclMatrix
                               jacobian, # vclVector  #form = c("loglik", "profileforBeta"),
                               minustwotimes=FALSE){
  
  
  ssqResidual <- as.matrix(ssqResidual)
  detVar <- as.vector(detVar)
  detVar <- matrix(rep(detVar, Ndata), nrow=Nparam)
  jacobian <- as.vector(jacobian)
  jacobian <- do.call(rbind, replicate(Nparam, jacobian, simplify = FALSE))
  m <- length(stderror)
  LogLik_optimized = matrix(0, nrow=m, ncol=1)
  
  
  
  for (var in 1:m){
    All_min2loglik_forthisvar <- ssqResidual/(stderror[var]^2) + Nobs*log(stderror[var]^2) + detVar + jacobian + Nobs*log(2*pi) 
    LogLik_optimized[var,] <- -0.5*min(All_min2loglik_forthisvar)
  }

  result = cbind(stderror, LogLik_optimized)
  colnames(result) <- c("stderror",'LogL')
  MLE <- result[,'stderror'][which.max(result[,'LogL'])]
  # plot(result[,'stderror'], result[,'LogL'])
  maximum <- max(LogLik_optimized)
  breaks95 = maximum - qchisq(0.95,  df = 1)/2
  
  
  leftOfMax = result[,'stderror'] < MLE
  if(length(which(leftOfMax)) <=2 | length(which(!leftOfMax)) <=2){
    ci95 = c(NA,NA)
    print("Not enough data for CI calculation")
  }else{
    afLeft <- approxfun(result[,'LogL'][leftOfMax], result[,'stderror'][leftOfMax])   
    afRight <- approxfun(result[,'LogL'][!leftOfMax], result[,'stderror'][!leftOfMax])   
  
    ci95= c(afLeft(breaks95), afRight(breaks95))
  }
    
  
  Output = list(LogLik_optimized=LogLik_optimized,
                  ci95=ci95
                  )
     
  Output
  
}


