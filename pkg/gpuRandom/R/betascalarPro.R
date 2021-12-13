#' @title Estimate Log-likelihood for Gaussian random fields when Betas are given
#' @useDynLib gpuRandom


     makeSymm <- function(m) {
         m[upper.tri(m)] <- t(m)[upper.tri(m)]
         return(m)
     }
     
     
#' @export
  betascalarProfile <- function(Betas, #a m x 1 R vector  given by the user 
                                cilevel=0.95,
                                a,     # which beta_i?
                                Nobs,  # number of observations.
                                Ndata,
                                Nparam,
                                Ncov,
                                detVar, # vclVector
                                ssqY,   # vclMatrix
                                XVYXVX,   # vclMatrix
                                jacobian # vclVector  #form = c("loglik", "profileforBeta"),
                                ){ 
  
  m <- length(Betas)
  if(m < 5){
    warning("need more values for accurate estimate")
  }
  
  
  ssqY <- as.matrix(ssqY)
  detVar <- as.vector(detVar)
  detVar <- matrix(rep(detVar, Ndata), nrow=Nparam)
  jacobian <- as.vector(jacobian)
  jacobian <- do.call(rbind, replicate(Nparam, jacobian, simplify = FALSE))
  #dim(XVYXVX)
  Ucov <- Ncov-1
  XTVinvX <- XVYXVX[ , (ncol(XVYXVX)-Ncov+1):ncol(XVYXVX)]
  XVY <- XVYXVX[ , 1:Ndata]
  
  #dim(XTVinvX)
  ## make each symmetric
  for (i in 1:Nparam){
    XTVinvX[((i-1)*Ncov+1) : (i*Ncov), ] <- makeSymm(XTVinvX[((i-1)*Ncov+1) : (i*Ncov), ])
    #mat <- XTVinvX[((i-1)*Ncov+1) : (i*Ncov), ]
    #XTVinvX[((i-1)*Ncov+1) : (i*Ncov), ] <-as.matrix(forceSymmetric(mat, "L"))
    #XTVinvX[((i-1)*Ncov+1) : (i*Ncov), ][upper.tri(XTVinvX[((i-1)*Ncov+1) : (i*Ncov), ])] <- XTVinvX[((i-1)*Ncov+1) : (i*Ncov), ][lower.tri(XTVinvX[((i-1)*Ncov+1) : (i*Ncov), ])]
  }
  
  #  XTVinvX[1:10,]
  # matrix(XTVinvX[-selectedrows,-a], ncol=Ucov, byrow=TRUE)[1:5,]
  # matrix(XTVinvX[selectedrows, -a], nrow=Nparam, ncol=Ucov)[1:5,]
  selectedrows <- (seq_len(Nparam)-1) * Ncov + a
  XTVinvX_deleted <- matrix(XTVinvX[-selectedrows,-a], ncol=Ucov, byrow=TRUE)
  XTVinvX_a <- matrix(XTVinvX[selectedrows, -a], nrow=Nparam, ncol=Ucov)
  #XTVinvX_a <- matrix(XTVinvX[-selectedrows, a], nrow=Nparam*Ucov, ncol=1)[1:10,]
  XVY_deleted <- XVY[-selectedrows, ]
  X_aVY <- XVY[selectedrows, ]
  X_aVX_a <- XTVinvX[selectedrows, a]
  #dim(XVY_deleted)
  
  
    partA = matrix(0, nrow=Nparam, ncol=Ndata)
    partB = matrix(0, nrow=Nparam, ncol=Ndata)
    partC = matrix(0, nrow=Nparam, ncol=Ndata)
    partD = matrix(0, nrow=Nparam, ncol=Ndata)
    partE = matrix(0, nrow=Nparam, ncol=Ndata)
    minus2LogLik_optimized = matrix(0, nrow=m, ncol=1)
    
    
  for (i in 1:Nparam){
    interval <- c(((i-1)*Ucov+1) : (i*Ucov))
    temp <- solve(XTVinvX_deleted[interval, ]) 
    
    # part (A) have 2 data sets
    partA[i,] = rowSums(XVY_deleted[interval,] %*% temp) * XVY_deleted[interval,]
    # part (B) have 2 data sets.   has beta
    partB[i,] = -2*XVY_deleted[interval,] %*% temp %*% XTVinvX_a[i,]
    # part (C) no data sets.  has beta
    partC[i,] = XTVinvX_a[i, ] %*% temp %*% XTVinvX_a[i, ]
    #partC[i,] = XTVinvX_a[interval, ] %*% temp %*% XTVinvX_a[interval, ]
    # part (D) have 2 data sets.    has beta
    partD[i,] = -2* X_aVY[i, ]
    # part (E)
    partE[i,] = X_aVX_a[i]
  }
    
    
    for (bet in 1:m){
    ssqResidual <- ssqY + Betas[bet] *partD + Betas[bet]^2 *partE - (partA + Betas[bet]* partB + Betas[bet]^2 * partC)
    All_min2loglik_forthisbeta <- Nobs*log(ssqResidual/Nobs) + detVar + Nobs + Nobs*log(2*pi) + jacobian
    minus2LogLik_optimized[bet,] = min(All_min2loglik_forthisbeta)
    }
    
    ############### ci ###########################################
    lower = min(Betas)
    upper = max(Betas)
    LogLik <- -0.5*minus2LogLik_optimized
    #f1 <- splinefun(Betas, LogLik, method = "fmm")
    f1 <- approxfun(Betas, LogLik)
    #plot(Betas,LogLik)
    #curve(f1(x), add = TRUE, col = 2, n = 1001)
    
    result <- optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)
    MLE <- result$maximum
    breaks <- result$objective - qchisq(cilevel,  df = 1)/2
    #abline(h=breaks)
    #f2 <- splinefun(Betas, LogLik-breaks, method = "fmm")
    f2 <- approxfun(Betas, LogLik-breaks)
    #plot(Betas,LogLik-breaks)
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
      warning("need wider range of beta to search ci's")
      ci <- c(NA, NA)
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
  
  
  
  
# likLgm_betascalarProfile <- function(Betas, #a m x 1 R vector  given by the user 
#                                      a,     # which beta_i?
#                                      Nobs,  # number of observations.
#                                      Ndata,
#                                      Nparam,
#                                      Ncov,
#                                      detVar, # vclVector
#                                      ssqY,   # vclMatrix
#                                      XVYXVX,   # vclMatrix
#                                      jacobian, # vclVector  #form = c("loglik", "profileforBeta"),
#                                      minustwotimes=TRUE){
#   
#   
#   ssqY <- as.matrix(ssqY)
#   detVar <- as.vector(detVar)
#   detVar <- matrix(rep(detVar, Ndata), nrow=Nparam)
#   jacobian <- as.vector(jacobian)
#   jacobian <- do.call(rbind, replicate(Nparam, jacobian, simplify = FALSE))
#   #dim(XVYXVX)
#   Ucov <- Ncov-1
#   m <- length(Betas)
#   
#   XTVinvX <- XVYXVX[ , (ncol(XVYXVX)-Ncov+1):ncol(XVYXVX)]
#   XVY <- XVYXVX[ , 1:Ndata]
#   
#   ## make each symmetric
#   # for (i in 1:Nparam){
#   #   #mat <- XTVinvX[((i-1)*Ncov+1) : (i*Ncov), ]
#   #   XTVinvX[((i-1)*Ncov+1) : (i*Ncov), ][upper.tri(XTVinvX[((i-1)*Ncov+1) : (i*Ncov), ])] <- XTVinvX[((i-1)*Ncov+1) : (i*Ncov), ][lower.tri(XTVinvX[((i-1)*Ncov+1) : (i*Ncov), ])]
#   # }
#   
#   
#   selectedrows <- seq_len(Nparam)*a
#   XTVinvX_deleted <- matrix(XTVinvX[-selectedrows,-a],ncol=Ucov, byrow=TRUE)
#   XTVinvX_a <- matrix(XTVinvX[selectedrows, -a], nrow=Nparam, ncol=Ucov)
#   #XTVinvX_a <- matrix(XTVinvX[-selectedrows, a], nrow=Nparam*Ucov, ncol=1)[1:10,]
#   XVY_deleted <- XVY[-selectedrows, ]
#   X_aVY <- XVY[selectedrows, ]
#   X_aVX_a <- XTVinvX[selectedrows, a]
#   #dim(XVY_deleted)
#   
#   
#   partA = matrix(0, nrow=Nparam, ncol=Ndata)
#   partB = matrix(0, nrow=Nparam, ncol=Ndata)
#   partC = matrix(0, nrow=Nparam, ncol=Ndata)
#   partD = matrix(0, nrow=Nparam, ncol=Ndata)
#   partE = matrix(0, nrow=Nparam, ncol=Ndata)
#   LogLik_optimized = matrix(0, nrow=m, ncol=1)
#   
#   
#   
#   
#   i_vec <- seq_len(Nparam)
#   computed <- parallel::mcmapply(i_vec,  
#                                  FUN = function(i) {
#                                    interval <- c(((i-1)*Ucov+1) : (i*Ucov))
#                                    solve(XTVinvX_deleted[interval, ]) 
#                                  }
#   )
#   
#   computed = matrix(computed, ncol = Nparam)
#   
#   for (i in 1:Nparam){
#     interval <- c(((i-1)*Ucov+1) : (i*Ucov))
#     temp <- matrix(computed[,i], ncol = Ucov)
#     
#     # part (A) have 2 data sets
#     partA[i,] = rowSums(XVY_deleted[interval,] %*% temp) * XVY_deleted[interval,]
#     # part (B) have 2 data sets.   has beta
#     partB[i,] = -2*XVY_deleted[interval,] %*% temp %*% XTVinvX_a[i,]
#     # part (C) no data sets.  has beta
#     partC[i,] = XTVinvX_a[i, ] %*% temp %*% XTVinvX_a[i, ]
#     #partC[i,] = XTVinvX_a[interval, ] %*% temp %*% XTVinvX_a[interval, ]
#     # part (D) have 2 data sets.    has beta
#     partD[i,] = -2* X_aVY[i, ]
#     # part (E)
#     partE[i,] = X_aVX_a[i]
#   }
#   
#   
#   for (bet in 1:m){
#     ssqResidual <- ssqY + Betas[bet] *partD + Betas[bet]^2 *partE - (partA + Betas[bet]* partB + Betas[bet]^2 * partC)
#     All_min2loglik_forthisbeta <- Nobs*log(ssqResidual/Nobs) + detVar + Nobs + Nobs*log(2*pi) + jacobian
#     LogLik_optimized[bet,] = min(All_min2loglik_forthisbeta)
#   }
#   
#   
#   LogLik_optimized
#   
# }  
