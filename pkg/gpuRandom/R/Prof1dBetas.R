#' @title Estimate Log-likelihood for Gaussian random fields when Betas are given
#' @useDynLib gpuRandom


     makeSymm <- function(m) {
         m[upper.tri(m)] <- t(m)[upper.tri(m)]
         return(m)
     }
  
#' @export

        Prof1dBetas <- function(Betas, #a m x p R matrix  given by the user 
                                cilevel,
                                Nobs,  # number of observations.
                                Ndata,
                                Nparam,
                                Ncov,
                                detVar, # 
                                ssqY,   # 
                                XVYXVX,   # 
                                jacobian){ 
  
  m <- nrow(Betas)
  # if(m < 5){
  #   warning("need more values for accurate estimate")
  # }
  

  detVar <- matrix(rep(detVar, Ndata), nrow=Nparam)
  jacobian <- do.call(rbind, replicate(Nparam, jacobian, simplify = FALSE))
  #dim(XVYXVX)
  Ucov <- Ncov-1
  XTVinvX <- XVYXVX[ , (ncol(XVYXVX)-Ncov+1):ncol(XVYXVX)]
  XVY <- matrix(XVYXVX[ , 1:Ndata], ncol=Ndata)
  
  #dim(XTVinvX)
  ## make each symmetric
  for (i in 1:Nparam){
    XTVinvX[((i-1)*Ncov+1) : (i*Ncov), ] <- makeSymm(XTVinvX[((i-1)*Ncov+1) : (i*Ncov), ])
  }
  
  LogLik_optimized = matrix(0, nrow=m, ncol=ncol(Betas))
  breaks <- rep(0, ncol(Betas))
  Table <- matrix(NA, nrow=ncol(Betas), ncol=4)
  colnames(Table) <-  c("MLE", "maximum", paste(c('lower', 'upper'), cilevel*100, 'ci', sep = ''))
  
  
  for (a in 1:ncol(Betas)){
  BetaSlice <- Betas[,a]
  selectedrows <- (seq_len(Nparam)-1) * Ncov + a
  XTVinvX_deleted <- matrix(XTVinvX[-selectedrows,-a], ncol=Ucov)
  XTVinvX_a <- matrix(XTVinvX[selectedrows, -a], nrow=Nparam, ncol=Ucov)
  XVY_deleted <- matrix(XVY[-selectedrows, ],ncol=Ndata)
  X_aVY <- matrix(XVY[selectedrows, ], ncol=Ndata)
  X_aVX_a <- XTVinvX[selectedrows, a]
  
  partA = matrix(0, nrow=Nparam, ncol=Ndata)
  partB = matrix(0, nrow=Nparam, ncol=Ndata)
  partC = matrix(0, nrow=Nparam, ncol=Ndata)
  partD = matrix(0, nrow=Nparam, ncol=Ndata)
  partE = matrix(0, nrow=Nparam, ncol=Ndata)
  
  if(Ndata == 1){
    for (i in 1:Nparam){
      interval <- c(((i-1)*Ucov+1) : (i*Ucov))
      eigenH = eigen(XTVinvX_deleted[interval, ])
      eig = list(values=1/eigenH$values, vectors=eigenH$vectors)
      temp <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
      #temp <- solve(XTVinvX_deleted[interval, ]) 
      #temp %*% XTVinvX_deleted[interval, ]
      # part (A) have 2 data sets
      partA[i,] = XVY_deleted[interval,] %*% temp %*% XVY_deleted[interval,]
      # part (B) have 2 data sets.   has beta
      partB[i,] = -2*XVY_deleted[interval,] %*% temp %*% XTVinvX_a[i,]
      # part (C) no data sets.  has beta
      partC[i,] = XTVinvX_a[i, ] %*% temp %*% XTVinvX_a[i, ]
      #partC[i,] = XTVinvX_a[interval, ] %*% temp %*% XTVinvX_a[interval, ]
      # part (D) have 2 data sets.    has beta
      partD[i,] = -2* X_aVY[i, ]
      # part (E)
      partE[i,] = X_aVX_a[i]
      #print(i)
    }
  }else{
    for (i in 1:Nparam){
      interval <- c(((i-1)*Ucov+1) : (i*Ucov))
      temp <- solve(XTVinvX_deleted[interval, ]) 
      
      if(Ucov==1){
        #diag(XVY_deleted[interval,] %*% temp  %*% XVY_deleted[interval,])
        partA[i,] = rowSums(XVY_deleted[interval,] %*% temp * XVY_deleted[interval,])
        # part (B) have 2 data sets.   has beta
        partB[i,] = -2*XVY_deleted[interval,] %*% temp %*% XTVinvX_a[i,]
        # part (C) no data sets.  has beta
        partC[i,] = XTVinvX_a[i, ] %*% temp %*% XTVinvX_a[i, ]
        #partC[i,] = XTVinvX_a[interval, ] %*% temp %*% XTVinvX_a[interval, ]
        # part (D) have 2 data sets.    has beta
        partD[i,] = -2* X_aVY[i, ]
        # part (E)
        partE[i,] = X_aVX_a[i]
      }else{
      #diag(t(XVY_deleted[interval,]) %*% temp %*% XVY_deleted[interval,])
      # part (A) have 2 data sets
      partA[i,] = rowSums(t(XVY_deleted[interval,]) %*% temp * t(XVY_deleted[interval,]))
      # part (B) have 2 data sets.   has beta
      partB[i,] = -2*t(XVY_deleted[interval,]) %*% temp %*% XTVinvX_a[i,]
      # part (C) no data sets.  has beta
      partC[i,] = XTVinvX_a[i, ] %*% temp %*% XTVinvX_a[i, ]
      #partC[i,] = XTVinvX_a[interval, ] %*% temp %*% XTVinvX_a[interval, ]
      # part (D) have 2 data sets.    has beta
      partD[i,] = -2* X_aVY[i, ]
      # part (E)
      partE[i,] = X_aVX_a[i]
      }
    }
  }
  
  # loglikAll = array(NA, c(m,Nparam,Ndata))
  # 
  # for (bet in 1:m){
  #     ssqResidual <- ssqY + Betas[bet] *partD + Betas[bet]^2 *partE - (partA + Betas[bet]* partB + Betas[bet]^2 * partC)
  #     loglik_forthisbeta <- (-0.5)*(Nobs*log(ssqResidual/Nobs) + detVar + Nobs + Nobs*log(2*pi) + jacobian)
  #     
  #     loglikAll[bet,,] = loglik_forthisbeta
  #     
  #     #LogLik_optimized[bet,] = max(loglik_forthisbeta)
  # }
  
#   plot(Betas, loglikAll[,1,12], col='red', lwd=2, type='l', ylim = c(-50,0) + max(loglikAll))
#   abline(v=6.04, col='blue')
#   
# bestParam =   which.max(apply(loglikAll, c(1,3), max))
# bestBocxox = which.max(apply(loglikAll[,bestParam,], 2, max))  
# lines(Betas, loglikAll[,bestParam, bestBocxox], col='blue')

  for (bet in 1:m){
    ssqResidual <- ssqY + BetaSlice[bet] *partD + BetaSlice[bet]^2 *partE - (partA + BetaSlice[bet]* partB + BetaSlice[bet]^2 * partC)
    loglik_forthisbeta <- (-0.5)*(Nobs*log(ssqResidual/Nobs) + detVar + Nobs + Nobs*log(2*pi) + jacobian)
    LogLik_optimized[bet,a] = max(loglik_forthisbeta)
  }
  

  
  # Sconfigs = c(1, seq(4*608, 608*10, len=50))
  # xx <- matrix(0, nrow=m, ncol=length(Sconfigs))
  # par(mfrow = c(1,1), mar = c(3,3, 0, 0))
  # for (bet in 1:m){
  #   ssqResidual <- ssqY + BetaSlice[bet] *partD + BetaSlice[bet]^2 *partE - (partA + BetaSlice[bet]* partB + BetaSlice[bet]^2 * partC)
  #   loglik_forthisbeta <- (Nobs*log(ssqResidual/Nobs) + detVar + Nobs + Nobs*log(2*pi) + jacobian)*(-0.5)
  #   loglik_forthisbeta <- loglik_forthisbeta[Sconfigs]
  #   
  #   xx[bet,] <- loglik_forthisbeta #[order(loglik_forthisbeta,decreasing = TRUE)]#[1:100]]
  # }
  # which (loglik_forthisbeta == min(loglik_forthisbeta))
  
#   matplot(BetaSlice, xx, type='l', xlab='intercept', lty=1, col=c("#FF0000", rep("#00000070", length(Sconfigs)-1)),
#           ylim = max(xx) + c(-6, 0))#, ylim = max(loglik_forthisbeta) + c(-1, 0))
#   lines(BetaSlice, xx[,1], col='red')
#        #lines(BetaSlice, loglik_forthisbeta[])
#   loglik_forthisbeta[order(loglik_forthisbeta,decreasing = TRUE)[1:100]]
  #  aaa<-matrix(loglikAll[,,12], nrow=m, ncol=Nparam)
  #  matplot(BetaSlice, aaa, type='l', xlab='intercept', lty=1, ylim=c(-395,-377),   col='#00000040')
  #  abline(v=c(simRes$summary['(Intercept)',c('estimate','ci0.1','ci0.9')]))
  #  lines(BetaSlice, LogLik_optimized, col='red')
  #  plot(BetaSlice, LogLik_optimized, col='red')
  #  dim(aaa)
  # 
  # aaa[,1:3]
  
  
    ############### ci ###########################################
    lower = min(BetaSlice)
    upper = max(BetaSlice)
    LogLik <- LogLik_optimized[,a]
    #f1 <- splinefun(Betas, LogLik, method = "fmm")
    breaks[a] <- max(LogLik) - qchisq(cilevel,  df = 1)/2
    f1 <- approxfun(BetaSlice, LogLik-breaks[a])
    #plot(BetaSlice,LogLik-breaks[a], cex=0.2, ylim=c(0.45,0.85), xlim=c(5,7.5))
    #curve(f1(x), add = TRUE, col = 2, n = 1001)
    
    result <- optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)
    MLE <- result$maximum
    
    plot(BetaSlice, LogLik-breaks[a],  ylim = max(LogLik-breaks[a]) + c(-3, 0.2), xlim = range(BetaSlice[max(LogLik-breaks[a]) - LogLik+breaks[a] < 3]), cex=0.2, xlab=paste('beta',a))
    abline(h=0, lty = 2, col=2)
    #abline(h=breaks[a], lty = 2)
    #f2 <- splinefun(Betas, LogLik-breaks[a], method = "fmm")
    #f2 <- approxfun(Betas, LogLik-breaks[a])
    #plot(Betas,LogLik-breaks[a])
    #curve(f2(x), add = TRUE, col = 2, n = 1001)
    ci <- rootSolve::uniroot.all(f1, lower = lower, upper = upper)
    abline(v=c(MLE,ci), lty=2)
    if(length(ci)==1){
      if( ci > MLE){
      ci <- c(lower, ci)
      }else{
      ci <- c(ci, upper)}
    }
    
    if(length(ci)==0 | length(ci)>2){
      warning("invalid ci returned")
      ci <- c(NA, NA)
    }

    ############### output #####################################
    Table[a,] <- c(MLE, result$objective+breaks[a], ci)
    
  }    
    Output <- list(estimates = Table,
                   LogLik = LogLik_optimized,
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
#     loglik_forthisbeta <- Nobs*log(ssqResidual/Nobs) + detVar + Nobs + Nobs*log(2*pi) + jacobian
#     LogLik_optimized[bet,] = min(loglik_forthisbeta)
#   }
#   
#   
#   LogLik_optimized
#   
# }  
        
        
        
        
        
        
        
        
        
        
        # for (i in 1:Nparam){
        #   XTVinvX[((i-1)*Ncov+1) : (i*Ncov), ] <- makeSymm(XTVinvX[((i-1)*Ncov+1) : (i*Ncov), ])
        #   #mat <- XTVinvX[((i-1)*Ncov+1) : (i*Ncov), ]
        #   #XTVinvX[((i-1)*Ncov+1) : (i*Ncov), ] <-as.matrix(forceSymmetric(mat, "L"))
        #   #XTVinvX[((i-1)*Ncov+1) : (i*Ncov), ][upper.tri(XTVinvX[((i-1)*Ncov+1) : (i*Ncov), ])] <- XTVinvX[((i-1)*Ncov+1) : (i*Ncov), ][lower.tri(XTVinvX[((i-1)*Ncov+1) : (i*Ncov), ])]
        # }
        # 
        # #  XTVinvX[1:10,]
        # # matrix(XTVinvX[-selectedrows,-a], ncol=Ucov, byrow=TRUE)[1:5,]
        # # matrix(XTVinvX[selectedrows, -a], nrow=Nparam, ncol=Ucov)[1:5,]
        # selectedrows <- (seq_len(Nparam)-1) * Ncov + a
        # XTVinvX_deleted <- matrix(XTVinvX[-selectedrows,-a], ncol=Ucov, byrow=TRUE)
        # XTVinvX_a <- matrix(XTVinvX[selectedrows, -a], nrow=Nparam, ncol=Ucov)
        # #XTVinvX_a <- matrix(XTVinvX[-selectedrows, a], nrow=Nparam*Ucov, ncol=1)[1:10,]
        # XVY_deleted <- matrix(XVY[-selectedrows, ], nrow = Nparam*Ucov, byrow=TRUE)
        # X_aVY <- matrix(XVY[selectedrows, ], nrow=Nparam, byrow=TRUE)
        # X_aVX_a <- XTVinvX[selectedrows, a]
        # #dim(XVY_deleted)
        # 
        # 
        # partA = matrix(0, nrow=Nparam, ncol=Ndata)
        # partB = matrix(0, nrow=Nparam, ncol=Ndata)
        # partC = matrix(0, nrow=Nparam, ncol=Ndata)
        # partD = matrix(0, nrow=Nparam, ncol=Ndata)
        # partE = matrix(0, nrow=Nparam, ncol=Ndata)
        # LogLik_optimized = matrix(0, nrow=m, ncol=1)
        # 
        # 
        # for (i in 1:Nparam){
        #   interval <- c(((i-1)*Ucov+1) : (i*Ucov))
        #   temp <- solve(XTVinvX_deleted[interval, ]) 
        #   
        #   # part (A) have 2 data sets
        #   partA[i,] = rowSums(XVY_deleted[interval,] %*% temp) * XVY_deleted[interval,]
        #   # part (B) have 2 data sets.   has beta
        #   partB[i,] = -2*XVY_deleted[interval,] %*% temp %*% XTVinvX_a[i,]
        #   # part (C) no data sets.  has beta
        #   partC[i,] = XTVinvX_a[i, ] %*% temp %*% XTVinvX_a[i, ]
        #   #partC[i,] = XTVinvX_a[interval, ] %*% temp %*% XTVinvX_a[interval, ]
        #   # part (D) have 2 data sets.    has beta
        #   partD[i,] = -2* X_aVY[i, ]
        #   # part (E)
        #   partE[i,] = X_aVX_a[i]
        # }
        # 
        # 
        # for (bet in 1:m){
        #   ssqResidual <- ssqY + Betas[bet] *partD + Betas[bet]^2 *partE - (partA + Betas[bet]* partB + Betas[bet]^2 * partC)
        #   loglik_forthisbeta <- Nobs*log(ssqResidual/Nobs) + detVar + Nobs + Nobs*log(2*pi) + jacobian
        #   LogLik_optimized[bet,] = min(loglik_forthisbeta)
        # }
        