#' @title Estimate Log-likelihood for Gaussian random fields when Betas are given
#' @useDynLib gpuRandom
#' @export



likLgm_betascalarProfile <- function(Betas, #a m x 1 R vector  given by the user 
                                     a,     # which beta_i?
                                     Nobs,  # number of observations.
                                     Ndata,
                                     Nparam,
                                     Ncov,
                                     detVar, # vclVector
                                     ssqY,   # vclMatrix
                                     XVYXVX,   # vclMatrix
                                     jacobian, # vclVector  #form = c("loglik", "profileforBeta"),
                                     minustwotimes=TRUE){
  

  ssqY <- as.matrix(ssqY)
  detVar <- as.vector(detVar)
  detVar <- matrix(rep(detVar, Ndata), nrow=Nparam)
  jacobian <- as.vector(jacobian)
  jacobian <- do.call(rbind, replicate(Nparam, jacobian, simplify = FALSE))
  #dim(XVYXVX)
  Ucov <- Ncov-1
  m <- length(Betas)
  
  XTVinvX <- XVYXVX[ , (ncol(XVYXVX)-Ncov+1):ncol(XVYXVX)]
  XVY <- XVYXVX[ , 1:Ndata]
  
  ## make each symmetric
  for (i in 1:Nparam){
    #mat <- XTVinvX[((i-1)*Ncov+1) : (i*Ncov), ]
    XTVinvX[((i-1)*Ncov+1) : (i*Ncov), ][upper.tri(XTVinvX[((i-1)*Ncov+1) : (i*Ncov), ])] <- XTVinvX[((i-1)*Ncov+1) : (i*Ncov), ][lower.tri(XTVinvX[((i-1)*Ncov+1) : (i*Ncov), ])]
  }
  
  
  selectedrows <- seq_len(Nparam)*a
  XTVinvX_deleted <- matrix(XTVinvX[-selectedrows,-a],ncol=Ucov, byrow=TRUE)
  XTVinvX_a <- matrix(XTVinvX[-selectedrows, a], nrow=Ucov*Nparam, ncol=1)
  XVY_deleted <- XVY[-selectedrows, ]
  X_aVY <- XVY[selectedrows, ]
  X_aVX_a <- XTVinvX[selectedrows, a]
  #dim(XVY_deleted)
  
  
    partA = matrix(0, nrow=Nparam, ncol=Ndata)
    partB = matrix(0, nrow=Nparam, ncol=Ndata)
    partC = matrix(0, nrow=Nparam, ncol=Ndata)
    partD = matrix(0, nrow=Nparam, ncol=Ndata)
    partE = matrix(0, nrow=Nparam, ncol=Ndata)
    LogLik_optimized = matrix(0, nrow=m, ncol=1)
    
    
  for (i in 1:Nparam){
    interval <- c(((i-1)*Ucov+1) : (i*Ucov))
    temp <- solve(XTVinvX_deleted[interval, ]) 
    
    # part (A) have 2 data sets
    partA[i,] = rowSums(XVY_deleted[interval,] %*% temp) * XVY_deleted[interval,]
    # part (B) have 2 data sets.   has beta
    partB[i,] = -2*XVY_deleted[interval,] %*% temp %*% XTVinvX_a[interval,]
    # part (C) no data sets.  has beta
    partC[i,] = XTVinvX_a[interval, ] %*% temp %*% XTVinvX_a[interval, ]
    # part (D) have 2 data sets.    has beta
    partD[i,] = -2* X_aVY[interval, ]
    # part (E)
    partE[i,] = X_aVX_a[i]
  }
    
    
    for (bet in 1:m){
    ssqResidual <- ssqY + Betas[bet] *partD + Betas[bet]^2 *partE - (partA + Betas[bet]* partB + Betas[bet]^2 * partC)
    All_min2loglik_forthisbeta <- Nobs*log(ssqResidual/Nobs) + detVar + Nobs + Nobs*log(2*pi) + jacobian
    LogLik_optimized[bet,] = min(All_min2loglik_forthisbeta)
    }
  
  
    LogLik_optimized
  
}  
  
  
  
  
  