#' @title Estimate Log-likelihood for Gaussian random fields when Betas are given
#' @import data.table
#' @useDynLib gpuRandom
#' @export



 betavectorProfile <- function(Betas, #a p x 1 R matrix  given by the user 
                               prof2list,
                               Nobs,  # number of observations.
                               Nparam,
                               boxcox,
                               detVar, # vclVector
                               ssqY,   # vclMatrix
                               XVYXVX,   # vclMatrix
                               jacobian, # vclVector  #form = c("loglik", "profileforBeta"),
                               cilevel){
  
   # prof2list = list(Betas1 <- Betas1,
   #                  Betas2 <- Betas2)
   # 
   # Betas = do.call(expand.grid, prof2list)  
  
  if(!is.matrix(Betas))
  Betas<-as.matrix(Betas)
  
  #boxcox = c(1, 0, setdiff(boxcox, c(0,1)))
  Ncov = ncol(Betas)
  Ndata = length(boxcox)
  m = nrow(Betas)
  #Nparam = nrow(params)
  
  
  
  detVar <- matrix(rep(detVar, Ndata), nrow=Nparam)
  jacobian <- do.call(rbind, replicate(Nparam, jacobian, simplify = FALSE))
  XTVinvX <- XVYXVX[ , (ncol(XVYXVX)-Ncov+1):ncol(XVYXVX)]
  bTDinva <- as.matrix(XVYXVX[ , 1:Ndata], nrow=Ncov, Ncol=Ndata)
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
      midItem[i,] <- Betas[beta, ] %*% as.matrix(bTDinva[((i-1)*Ncov+1) : (i*Ncov), ], nrow=Ncov)
      
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
  
  LogLikBetas = -0.5*likForBetas
  result = data.frame(cbind(Betas,LogLikBetas))
  colnames(result) <- c('Betas1','Betas2',"LogLik")
  
  maximum <- max(LogLikBetas)
  breaks = maximum - qchisq(cilevel,  df = Ncov)/2
 
  lMatrix = matrix(result[,'LogLik'], length(prof2list[[1]]), length(prof2list[[2]]))
  contour(prof2list[[1]], prof2list[[2]], lMatrix,
          col = par("fg"), lty = par("lty"), lwd = par("lwd"), #levels = c(breaks), 
          add = FALSE, xlab = "Beta1", ylab = "Beta2")
  
  # filled.contour(prof2list[[1]], prof2list[[2]], lMatrix, levels = c(breaks-4, breaks-3, 
  #                                                                    breaks-2, breaks-1, breaks),
  #                color = function(n) hcl.colors(n, "ag_Sunset"),
  #                plot.title={
  #                  title(xlab = "Intercept", ylab = "Elevation",main = "contour plot")
  #                  #abline(h=trueParam['nugget'], col='red')
  #                  #abline(v=trueParam['range']/100, col = "red")
  #                  })
  

  Sprob = c(0,0.2, 0.5, 0.8, 0.9, 0.95, 0.999)
  likCol = mapmisc::colourScale(drop(result$LogLik), breaks=max(lMatrix, na.rm=TRUE)- qchisq(Sprob, df=2), style='fixed', col='Spectral')
  points(x = result[,1], y = result[,2], pch=20, col=likCol$plot, cex=0.8)
  points(result[which.max(result$LogLik),],pch=20)
  mapmisc::legendBreaks('bottomright', breaks = rev(Sprob), col=likCol$col)

  # resultbeta1 = data.table::as.data.table(cbind(-0.5*LogLikBetas, Betas[,1]))
  # colnames(resultbeta1) <- c("LogLik", "Beta1")
  # profileLogLik_beta1 <- resultbeta1[, .(profile=max(.SD)), by=Beta1]
  # 
  # 
  # resultintercept = data.table::as.data.table(cbind(-0.5*LogLikBetas, Betas[,1]))
  # colnames(resultintercept) <- c("LogLik", "intercept")
  # profileLogLik_intercept <- resultintercept[, .(profile=max(.SD)), by=intercept]
  
  
  #plot(profileLogLik_beta1$Beta1,profileLogLik_beta1$profile)
  #plot(profileLogLik_intercept$intercept,profileLogLik_intercept$profile)
  
  Theoutput <- list(dataforplot=result,
                    breaks =breaks,
                    ssqForBetas = ssqForBetas)
  
  
  Theoutput
  
  
}
  
  
  




















  
  