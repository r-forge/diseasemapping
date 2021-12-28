#' @title Estimate profile Log-likelihood for covariance parameters and lambda
#' @import data.table
#' @importFrom rootSolve uniroot.all
#' @useDynLib gpuRandom
#' @export



# betahat and sigmahat 
# profile log-likelihood for each covariance parameters + lambda
 likfitLgmCov2d <- function(data,
                         formula, 
                         coordinates,
                         params, # CPU matrix for now, users need to provide proper parameters given their specific need
                         boxcox,  # boxcox is always estimated
                         cilevel=0.95,  # decimal
                         type = c("float", "double")[1+gpuInfo()$double_support],
                         reml=FALSE, 
                         NparamPerIter,
                         Nglobal,
                         Nlocal,
                         NlocalCache,
                         verbose=FALSE){
  
  # will always be c(1,0,.......)
  boxcox = c(1, 0, setdiff(boxcox, c(0,1)))
  

  
  
  # get rid of NAs in data
  data = data.frame(data)
  theNA = apply(  data[,all.vars(formula),drop=FALSE],
                  1, 
                  function(qq) any(is.na(qq))
  )
  noNA = !theNA
  
  
  ############## the X matrix #################################
  covariates = model.matrix(formula, data[noNA,])
  observations = all.vars(formula)[1]
  observations = data.matrix(data[noNA, observations, drop=FALSE])
  
  Nobs = nrow(data)
  Ncov = ncol(covariates)
  Ndata = length(boxcox)
  Nparam = nrow(params)
  
  # whole data set including columns for transformed y's
  yx = vclMatrix(cbind(observations, matrix(0, Nobs, Ndata-1), covariates), 
                 type=type)
  
  
  # coordinates
  coordsGpu = vclMatrix(coordinates, type=type)  
  # box-cox
  boxcoxGpu = vclVector(boxcox, type=type)
  
  # prepare params, make sure variance=1 in params
  params[,"variance"]=1 
  params0 = geostatsp::fillParam(params)
  paramsGpu = vclMatrix(cbind(params0, matrix(0, nrow(params0), 22-ncol(params0))),type=type)
  
  varMat = vclMatrix(0, Nobs*NparamPerIter, Nobs, type=type)
  cholDiagMat = vclMatrix(0, NparamPerIter, Nobs, type=type)
  ssqY <- vclMatrix(0, Nparam, Ndata, type=type)
  XVYXVX = vclMatrix(0, Nparam * Ncov, ncol(yx), type=type)
  ssqBetahat <- vclMatrix(0, Nparam, Ndata, type=type)
  detVar = vclVector(0, Nparam,type=type)
  detReml = vclVector(0, Nparam, type=type)
  jacobian = vclVector(0, Ndata, type=type)   
  ssqYX = vclMatrix(0, ncol(yx) * NparamPerIter, ncol(yx), type=type)
  ssqYXcopy = vclMatrix(0, ncol(yx) * NparamPerIter, ncol(yx), type=type)
  LinvYX = vclMatrix(0, Nobs * NparamPerIter, ncol(yx), type=type)
  QinvSsqYx = vclMatrix(0, NparamPerIter*Ncov, Ndata, type = type)
  cholXVXdiag = vclMatrix(0, NparamPerIter, Ncov, type=type)
  minusTwoLogLik <- vclMatrix(0, Nparam, Ndata, type=type)
  
  
  gpuRandom:::likfitGpu_BackendP(
    yx,        #1
    coordsGpu, #2
    paramsGpu, #3
    boxcoxGpu,   #4#betasGpu,  #5
    ssqY,     #5#aTDinvb_beta,
    XVYXVX,
    ssqBetahat, #7#ssqBeta,
    detVar,    #11
    detReml,   #12
    jacobian,  #13
    NparamPerIter,  #14
    as.integer(Nglobal),  #12
    as.integer(Nlocal),  #16
    NlocalCache,  #14
    verbose=verbose,  #15
    ssqYX, #
    ssqYXcopy,  #new
    LinvYX,  #18
    QinvSsqYx, 
    cholXVXdiag, #20
    varMat,        #21     Vbatch
    cholDiagMat)
  
  # resid^T V^(-1) resid, resid = Y - X betahat = ssqResidual
  ssqResidual <- ssqY - ssqBetahat
  # any(is.na(as.matrix(log(ssqResidual/Nobs))))
   # as.matrix(cholDiagMat)[]
   # any(is.nan(as.matrix(log(ssqResidual/Nobs))))
   # any(is.nan(as.vector(detVar)))
   # any(is.na(as.matrix(varMat)))
   # any(is.nan(as.matrix(varMat)))
   # any(is.na(as.matrix(ssqYX)))
   # any(is.na(as.matrix(ssqY)))
   # any(is.na(as.matrix(ssqYXcopy)))
   # any(is.na(as.matrix(ssqBetahat)))
   # any(is.na(as.vector(detReml)))
   # any(is.na(as.matrix(cholXVXdiag)))
   # any(is.na(as.matrix(ssqResidual)))
   # any(is.nan(as.matrix(ssqResidual/Nobs)))
   # as.matrix(ssqY)[which(is.na(as.vector(detVar))),]
  # which(is.na(as.vector(detVar)))
  # debug <- as.matrix(log(ssqResidual/Nobs))
  # any(is.na(debug))
  # which(is.na(debug),arr.ind = TRUE)
  # params0[which(is.na(debug),arr.ind = TRUE)[1],]
  
  # ssqResidual[which(is.na(debug),arr.ind = TRUE)[,1],]     # row 20803 col 1  -4.530054e+20
  # ssqY[which(is.na(debug),arr.ind = TRUE)[,1],]
  # ssqBetahat[which(is.na(debug),arr.ind = TRUE)[,1],]    # problem arises from here! 4.530057e+20
  # as.vector(detVar)[840]
  # params0[which(is.na(as.vector(detVar))),]
  # params[which(is.na(result1$ssqY),arr.ind = TRUE)[1],]
  # if(fixVariance >0){
  #   temp <- vclMatrix(0, Nparam, Ndata, type=type)
  #   variances <- vclVector(params0[,3],type=type)
  #   mat_vec_eledivideBackend(ssqResidual, variances, temp,  Nglobal)
  # }
  
  if(reml== FALSE){ 
    # if(fixVariance == FALSE){ # ml
    gpuRandom:::matrix_vector_sumBackend(Nobs*log(ssqResidual/Nobs),
                                         detVar,
                                         jacobian,  
                                         Nobs + Nobs*log(2*pi),
                                         minusTwoLogLik,
                                         Nglobal)
    # }else{  #fixVariance == TRUE
    #   matrix_vector_sumBackend(temp,
    #                            Nobs*log(variances)+detVar,
    #                            jacobian,
    #                            Nobs*log(2*pi),
    #                            minusTwoLogLik,
    #                            Nglobal)
    # }
  }else if(reml==TRUE){
    # if(fixVariance == FALSE){# remlpro
    # minusTwoLogLik= (n-p)*log(ssqResidual/(n-p)) + logD + logP + jacobian + n*log(2*pi) + n-p    
    gpuRandom:::matrix_vector_sumBackend((Nobs-Ncov)*log(ssqResidual/(Nobs-Ncov)),
                                         detVar+detReml,
                                         jacobian,  
                                         Nobs*log(2*pi)+Nobs-Ncov,
                                         minusTwoLogLik,
                                         Nglobal)
    # }else{
    #   matrix_vector_sumBackend(temp,
    #                            detVar+detReml+(Nobs-Ncov)*log(variances),
    #                            jacobian,
    #                            Nobs*log(2*pi),
    #                            minusTwoLogLik,
    #                            Nglobal)
    #   
    #   }
  }
  
  
  LogLikcpu <- as.matrix(-0.5*minusTwoLogLik)
  maximum <- max(LogLikcpu)
  breaks = maximum - qchisq(cilevel,  df = 2)/2
  
  # if(twoDbreak == TRUE){
  #   breaks2D = maximum - qchisq(cilevel,  df = 2)/2 
  # }
  
  
  
  
  ############## output matrix ####################
  # Table <- matrix(NA, nrow=length(union(paramToEstimate, 'boxcox'))+Ncov+1, ncol=3)
  # rownames(Table) <-  c("intercept", paste(c('betahat'), seq_len(Ncov-1),sep = ''), "sdSpatial", union(paramToEstimate, 'boxcox'))
  # colnames(Table) <-  c("estimate", paste(c('lower', 'upper'), cilevel*100, 'ci', sep = ''))


  
  # index <- which(LogLikcpu == max(LogLikcpu, na.rm = TRUE), arr.ind = TRUE)
  # 
  # 
  # ###############betahat#####################
  # Betahat <- matrix(0, nrow=Ncov, ncol=Ndata)
  # a<-c( ((index[1]-1)*Ncov+1) : (index[1]*Ncov) )
  # mat <- XVYXVX[a,((Ndata+1):ncol(yx))]
  # mat[upper.tri(mat)] <- mat[lower.tri(mat)]
  # Betahat <- solve(mat) %*% XVYXVX[a,index[2]]
  # 
  # Table[c("intercept", paste(c('betahat'), seq_len(Ncov-1),sep = '')),1] <- Betahat
  # 
  # 
  # 
  # 
  # #################sigma hat#########################
  # if(reml==FALSE)  {
  #   Table["sdSpatial",1] <- sqrt(ssqResidual[index[1],index[2]]/Nobs)
  # }else{         
  #   Table["sdSpatial",1] <- sqrt(ssqResidual[index[1],index[2]]/(Nobs - Ncov))
  # }
  
  

  Output <- list(LogLik=LogLikcpu,
                 minusTwoLogLik = minusTwoLogLik,
                 breaks = breaks,
                 Nobs = Nobs,
                 Ncov = Ncov,
                 Ndata = Ndata,
                 Nparam = Nparam,
                 ssqY=ssqY,     
                 ssqBetahat = ssqBetahat,
                 ssqResidual = ssqResidual,
                 detVar = detVar,   
                 detReml = detReml,   
                 jacobian = jacobian,
                 XVYXVX=XVYXVX)
  
  
  
  Output
  
  
  
}






# Betahat <- matrix(0, nrow=Ncov, ncol=Ndata)
# # a<-c((index[1]-1)*Ncov+1, index[1]*Ncov)
# # Betahat <- solve(XVYXVX[a,((Ndata+1):ncol(yx))]) %*% XVYXVX[a,index[2]]
# for (j in 1:Ndata){
#   index2 <- order(LogLikcpu[,j], decreasing = TRUE)[1]  #from small to large
#   a<-c((index2[1]-1)*Ncov+1, index2[1]*Ncov)
#   Betahat[,j] <- solve(XVYXVX[a,((Ndata+1):ncol(yx))]) %*% XVYXVX[a,j]
# }
# Table[paste(c('betahat'),seq_len(Ncov),sep = ''),1] <- Betahat[,index[2]]
# 
# 












