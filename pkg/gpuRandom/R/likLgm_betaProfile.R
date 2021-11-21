#' @title Estimate Log-likelihood for Gaussian random fields when Betas are given
#' @useDynLib gpuRandom
#' @export



likLgm_betaProfile <- function(Betas, #a m x p R matrix  given by the user 
                               Nobs,  # number of observations.
                               Nparam,
                               BoxCox,
                               detVar, # vclVector
                               ssqY,   # vclMatrix
                               XVYXVX,   # vclMatrix
                               jacobian, # vclVector  #form = c("loglik", "profileforBeta"),
                               minustwotimes=TRUE){
  
  
  if(!is.matrix(Betas))
  Betas<-as.matrix(Betas)
  
  BoxCox = c(1, 0, setdiff(BoxCox, c(0,1)))
  Ncov = ncol(Betas)
  Ndata = length(BoxCox)
  m = nrow(Betas)
  #Nparam = nrow(params)
  
  
  
  ssqY <- as.matrix(ssqY)
  detVar <- as.vector(detVar)
  detVar <- matrix(rep(detVar, Ndata), nrow=Nparam)
  jacobian <- as.vector(jacobian)
  jacobian <- do.call(rbind, replicate(Nparam, jacobian, simplify = FALSE))
  #XVYXVX <- as.matrix(XVYXVX)
  XTVinvX <- XVYXVX[ , (ncol(XVYXVX)-Ncov+1):ncol(XVYXVX)]
  # bTDinva
  bTDinva <- XVYXVX[ , 1:Ndata]
  midItem <- matrix(88, nrow=Nparam, ncol=Ndata)
  ssqBeta0 <- rep(0, length=Nparam)   # dosen't depend on lambda
  ssqForBetas <- matrix(0, nrow=m, ncol=1)
  likForBetas<- matrix(0, nrow=m, ncol=1)
  
  
  
  # one = ssqY - 2*beta^T * bTDinva + ssqBeta
  # n*log(one/n)+ logD +jacobian + n*log(2*pi) + n
  # takes about over 7 minutes for 126,000 params
  for(beta in 1:m){  
    
    for (i in 1:Nparam){
      # for each of the beta, calculate the beta^T * bTDinva, 1 x Ndata
      midItem[i,] <- Betas[beta, ] %*% bTDinva[((i-1)*Ncov+1) : (i*Ncov), ]
      
      # ssqBeta = beta^T * (b^T D^(-1) b) * beta
      mat <- XTVinvX[((i-1)*Ncov+1) : (i*Ncov), ]
      mat[upper.tri(mat)] <- mat[lower.tri(mat)]
      ssqBeta0[i] <- Betas[beta, ] %*% mat %*% (Betas[beta, ])
    }
    
    ssqBeta <- do.call(cbind, replicate(Ndata, ssqBeta0, simplify=FALSE))   # ssq same for different data sets
    one <- ssqY - 2*midItem + ssqBeta
    
    temp <- Nobs*log(one/Nobs) + detVar + Nobs + Nobs*log(2*pi) + jacobian
    likForBetas[beta,] = min(temp)
    ssqForBetas[beta,] <- ssqBeta0[which(temp == min(temp, na.rm = TRUE), arr.ind = TRUE)[1]]
    
  }
  
  Theoutput <- list(likForBetas=likForBetas,
                    ssqForBetas = ssqForBetas)
  
  
  Theoutput
  
  
}
  
  
  
  
  
  