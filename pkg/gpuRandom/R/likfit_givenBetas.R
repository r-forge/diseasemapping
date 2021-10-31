#' @title Estimate Log-likelihood for Gaussian random fields when Betas are given
#' @useDynLib gpuRandom
#' @export



likfit_givenBeta <- function(Betas, #a p x m matrix  given by the user 
                             #DatasetIndex=NA,     #an integer, which dataset are the betas for?,
                             variances,
                             Nobs,  # number of observations.
                             detVar,
                             ssqY,
                             XVYXVX,
                             jacobian,
                             form = c("loglik", "mlFixSigma"),
                             minustwotimes=TRUE){

  form = c(loglik=1, mlFixSigma=2)[form]

  Ncov = nrow(Betas)
  m = ncol(Betas)   # number of sets of Betas
  Nparam = length(detVar)                          
  Ndata = length(jacobian)

  
  ssqY <- as.matrix(ssqY)
  detVar <- as.vector(detVar)
  detVar <- matrix(rep(detVar, Ndata), nrow=Nparam)
  jacobian <- as.vector(jacobian)
  jacobian <- do.call(rbind, replicate(Nparam, jacobian, simplify = FALSE))
  XVYXVX <- as.matrix(XVYXVX)
  
  

  XTVinvX <- XVYXVX[ , (ncol(XVYXVX)-Ncov+1):ncol(XVYXVX)]
  ssqForBeta <- matrix(0, nrow=Nparam, ncol=m)
  #maximized over lambda
  minusTwoLogLikOverLambda <- matrix(0, nrow=Nparam, ncol=m)    #row is param, col is Betas
  index <- matrix(0, nrow=m, ncol=2) #row is Betas, col is row & col index 
  
  
  for(beta in 1:m){  
    
  # param_vec  <- rep(seq_len(Nparam), each = Ndata)
  # lambda_vec <- rep(seq_len(Ndata), times = Nparam)
  # 
  # midItem <- parallel::mcmapply(param_vec, lambda_vec, 
  #                     FUN = function(i, j) {
  #                       # This is where expensive operations should go
  #                       t(XVYXVX[((i-1)*Ncov+1) : (i*Ncov), j] ) %*% Betas[ ,beta] 
  #                     }
  #                     )
  # 
  # # computed was a vector, but we need to put it in the correct-size matrix
  # midItem <- matrix(midItem, ncol = Ndata)
    
    
    midItem <- matrix(-66, nrow=Nparam, ncol=Ndata)
    # to calculate aTDinvb * Beta    Nparam by m
     for(lambda in 1:Ndata){
     for(i in 1:Nparam){
       midItem[i,lambda] <- t(XVYXVX[((i-1)*Ncov+1) : (i*Ncov), lambda]) %*% Betas[ ,beta]    # to check
     }
     }
     


  # ssqBeta = beta^T * (b^T D^(-1) b) * beta
  # param_vec  <- seq_len(Nparam)
  # ssqBeta0 <- parallel::mcmapply(param_vec, 
  #                FUN = function(i) {
  #                  t(Betas[ ,beta]) %*% XTVinvX[((i-1)*Ncov+1) : (i*Ncov), ] %*% Betas[ ,beta]  
  #                }
  #                )
    
    
    ssqBeta0 <- seq(0, length=Nparam)   # dosen't depend on lambda
    # ssqBeta = beta^T * (b^T D^(-1) b) * beta
     for (i in 1:Nparam){
            ssqBeta0[i] <- t(Betas[ ,beta]) %*% XTVinvX[((i-1)*Ncov+1) : (i*Ncov), ] %*% Betas[ ,beta]
        }
    
  
  ssqForBeta[,beta] <- ssqBeta0
  
  ssqBeta <- do.call(cbind, replicate(Ndata, ssqBeta0, simplify=FALSE))   # ssq for this Beta only
  
  
  one <- ssqY - 2*midItem + ssqBeta
  
  

    if(form == 1) {
      if (is.null(variances))
        stop("variances must be given")
      else{
        variances <- do.call(cbind, replicate(Ndata, variances, simplify=FALSE)) 
        temp <- Nobs*log(variances) + detVar + one/variances + Nobs*log(2*pi) + jacobian
        index[beta,] <- which(temp == min(temp), arr.ind = TRUE)
        minusTwoLogLikOverLambda[, beta] <-  apply( Nobs*log(variances) + detVar + one/variances + Nobs*log(2*pi) + jacobian, 1, min)
      }
    }
  
  
    if(form == 2){
      temp <- Nobs*log(one/Nobs) + detVar + Nobs + Nobs*log(2*pi) + jacobian
      index[beta,] <- which(temp == min(temp, na.rm = TRUE), arr.ind = TRUE)
      minusTwoLogLikOverLambda[, beta] <-  apply( temp, 1, min, na.rm=TRUE)
    }
  
  
  } 
  
    likForBeta = apply( minusTwoLogLikOverLambda, 2, min )
  
    Theoutput <- list(likForBeta=likForBeta,
                      minusTwoLogLikOverLambda = minusTwoLogLikOverLambda, 
                      index = index,
                      ssqForBeta = ssqForBeta)
    
    
    Theoutput
  
}











#midItem <- matrix(-66, nrow=Nparam, ncol=Ndata)
#ssqBeta0 <- vector(-88, length=Nparam)    # dosen't depend on lambda


## to calculate aTDinvb * Beta    Nparam by m
# for(beta in 1:m){
# for(lambda in 1:Ndata){
#   bTDinva <- XVYXVX[ , lambda] 
# for(i in 1:Nparam){     
#   midItem[i,lambda] <- t(bTDinva[((i-1)*Ncov+1) : (i*Ncov)]) %*% Betas[ ,beta]    # to check
# }
# } 
# }

## to calculate ssqBeta
# ssqBeta = beta^T * (b^T D^(-1) b) * beta
# for (i in 1:Nparam){
#        ssqBeta0[i] <- t(Betas[ ,beta]) %*% XTVinvX[((i-1)*Ncov+1) : (i*Ncov), ] %*% Betas[ ,beta]
#    }










