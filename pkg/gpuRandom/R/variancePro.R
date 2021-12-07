#' @title Estimate Log-likelihood for Gaussian random fields when stderror are given
#' @import data.table
#' @useDynLib gpuRandom
#' @export


######## this function will be incorporated in the other function later i think

       variancePro <- function(stderror, #a vector  given by the user 
                               cilevel=0.95,
                               Nobs,  # number of observations.
                               Nparam,
                               Ndata,
                               detVar, # vclVector
                               ssqResidual,   # vclMatrix
                               jacobian # vclVector  #form = c("loglik", "profileforBeta"),
                               ){
  
  
  ssqResidual <- as.matrix(ssqResidual)
  detVar <- as.vector(detVar)
  detVar <- matrix(rep(detVar, Ndata), nrow=Nparam)
  jacobian <- as.vector(jacobian)
  jacobian <- do.call(rbind, replicate(Nparam, jacobian, simplify = FALSE))
  m <- length(stderror)
  LogLik = matrix(0, nrow=m, ncol=1)
  
  
  
  for (var in 1:m){
    All_min2loglik_forthisvar <- ssqResidual/(stderror[var]^2) + Nobs*log(stderror[var]^2) + detVar + jacobian + Nobs*log(2*pi) 
    LogLik[var,] <- -0.5*min(All_min2loglik_forthisvar)
  }

  
  lower = min(stderror)
  upper = max(stderror)
  f1 <- splinefun(stderror, LogLik, method = "monoH.FC")
  #plot(stderror,LogLik)
  #curve(f1(x), add = TRUE, col = 2, n = 1001)
  
  result <- optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)
  MLE <- result$maximum
  breaks <- result$objective - qchisq(cilevel,  df = 1)/2
  #abline(h=breaks)
  f2 <- splinefun(stderror, LogLik-breaks, method = "monoH.FC")
  #plot(stderror,LogLik-breaks)
  #curve(f2(x), add = TRUE, col = 2, n = 1001)
  #abline(h=0)
  ci <- rootSolve::uniroot.all(f2, lower = lower, upper = upper)
  
  if(length(ci)==1){
    if( ci > MLE){
      ci <- c(lower, ci)
    }else{
      ci <- c(ci, upper)}
  }
  if(length(ci)==0){
    warning("require a better param matrix")
    ci <- c(NA, NA)
  }
  if(length(ci)>2){
    warning("need wider range of beta to search ci's")
    ci <- ci[3:4]
  }
  
  ############### output #####################################
  Table <- matrix(NA, nrow=1, ncol=4)
  colnames(Table) <-  c("MLE", "maximum", paste(c('lower', 'upper'), cilevel*100, 'ci', sep = ''))
  Table[1,] <- c(MLE, result$objective, ci)
  
  
  Output <- list(estimates = Table,
                 LogLik = LogLik,
                 breaks = breaks)
  
  Output
  
}

       

       # result = cbind(stderror, LogLik_optimized)
       # colnames(result) <- c("stderror",'LogL')
       # MLE <- result[,'stderror'][which.max(result[,'LogL'])]
       # # plot(result[,'stderror'], result[,'LogL'])
       # maximum <- max(LogLik_optimized)
       # breaks95 = maximum - qchisq(0.95,  df = 1)/2
       # 
       # 
       # leftOfMax = result[,'stderror'] < MLE
       # if(length(which(leftOfMax)) <=2 | length(which(!leftOfMax)) <=2){
       #   ci95 = c(NA,NA)
       #   print("Not enough data for CI calculation")
       # }else{
       #   afLeft <- approxfun(result[,'LogL'][leftOfMax], result[,'stderror'][leftOfMax])   
       #   afRight <- approxfun(result[,'LogL'][!leftOfMax], result[,'stderror'][!leftOfMax])   
       #   
       #   ci95= c(afLeft(breaks95), afRight(breaks95))
       # }
       
       
       # Output = list(LogLik_optimized=LogLik_optimized,
       #               ci95=ci95
       # )
       # 
       # Output
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       