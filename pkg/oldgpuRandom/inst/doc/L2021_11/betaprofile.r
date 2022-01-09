


Betas <- seq(-6, 15, len=20)*1e-04

a=2
cilevel=0.9
Nobs = result1$Nobs
Ndata = result1$Ndata
Nparam = nrow(result1$param)
Ncov = result1$Ncov
detVar = result1$detVar
ssqY = result1$ssqY
XVYXVX = result1$XVYXVX
jacobian = result1$jacobian

result1$ssqResidual

gpuRandom::Prof1dBetas(Betas, cilevel=0.9,a=2,  Nobs,  
                        Ndata,
                        Nparam,
                        Ncov,
                        detVar, 
                        ssqY,   
                        XVYXVX,   
                        jacobian)


m <- length(Betas)
if(m < 5){
  warning("need more values for accurate estimate")
}


detVar <- matrix(rep(detVar, Ndata), nrow=Nparam)
jacobian <- do.call(rbind, replicate(Nparam, jacobian, simplify = FALSE))
#dim(XVYXVX)
Ucov <- Ncov-1
XTVinvX <- XVYXVX[ , (ncol(XVYXVX)-Ncov+1):ncol(XVYXVX)]
XVY <- XVYXVX[ , 1:Ndata]

#dim(XTVinvX)
## make each symmetric
for (i in 1:Nparam){
  XTVinvX[((i-1)*Ncov+1) : (i*Ncov), ] <- makeSymm(XTVinvX[((i-1)*Ncov+1) : (i*Ncov), ])
}

selectedrows <- (seq_len(Nparam)-1) * Ncov + a
XTVinvX_deleted <- matrix(XTVinvX[-selectedrows,-a], ncol=Ucov, byrow=TRUE)
XTVinvX_a <- matrix(XTVinvX[selectedrows, -a], nrow=Nparam, ncol=Ucov)
XVY_deleted <- XVY[-selectedrows, ]
X_aVY <- XVY[selectedrows, ]
X_aVX_a <- XTVinvX[selectedrows, a]

partA = matrix(0, nrow=Nparam, ncol=Ndata)
partB = matrix(0, nrow=Nparam, ncol=Ndata)
partC = matrix(0, nrow=Nparam, ncol=Ndata)
partD = matrix(0, nrow=Nparam, ncol=Ndata)
partE = matrix(0, nrow=Nparam, ncol=Ndata)

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
#abline(h=breaks, lty = 2)
#f2 <- splinefun(Betas, LogLik-breaks, method = "fmm")
f2 <- approxfun(Betas, LogLik-breaks)
ci <- rootSolve::uniroot.all(f2, lower = lower, upper = upper)
