#' @title Estimate Log-likelihood for Gaussian random fields when Betas are given
#' @useDynLib gpuRandom
#' @export



likfit_givenBeta <- function(Betas, #a p x m matrix  given by the user 
                             Nparam,
                             DatasetIndex=NA,     #an integer, which dataset are the betas for?,
                             Nobs,  # number of observations.
                             detVar,
                             ssqY,
                             XVYXVX,
                             minustwotimes=TRUE){
  
  
  Ncov = nrow(Betas)
  
  m = ncol(Betas)
  
  detVar <- as.vector(derVar)
  detVar <- matrix(rep(detVar, m), nrow=Nparam)
  aTDinva <- matrix(rep(ssqY[,DatasetIndex],m),nrow=Nparam) 
  XTVinvX <- XVYXVX[ , (ncol(XVYXVX)-Ncov+1):ncol(XVYXVX)]
  bTDinva <- XVYXVX[ , DatasetIndex]
  
  
  middleItem <- matrix(-66, nrow=Nparam, ncol=m)
  ssqBeta <- matrix(-88, nrow=Nparam, ncol=m)
  
  
  ## to calculate aTDinvb * Beta    Nparam by m
  for(i in 1:Nparam){     
    middleItem[i,] <- t(bTDinva[((i-1)*Ncov+1) : (i*Ncov)]) %*% Betas    # to check
  }
  
  
  ## to calculate ssqBeta
  for (i in 1:Nparam){
    XTVinvX_i <- XTVinvX[((i-1)*Ncov+1) : (i*Ncov), ]
    
    for(j in 1:m){
      ssqBeta[i,j] <- t(Betas[ ,j]) %*% XTVinvX_i %*% Betas[ ,j]
    }
  }
   
  
  one <- aTDinva - 2*middleItem + ssqBeta
  
  minusTwoLogLik <-  Nobs*log(one/Nobs) + detVar + Nobs + Nobs*log(2*pi)
  
  
}


