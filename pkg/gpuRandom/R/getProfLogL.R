#' @title gpu elements
#' @useDynLib gpuRandom
#' @export





  getProfLogL <- function(data,
                          formula, 
                          coordinates,
                          params,  # CPU matrix 
                          boxcox,  # boxcox is always estimated
                          type = c("float", "double")[1+gpuInfo()$double_support],
                          NparamPerIter,
                          gpuElementsOnly = FALSE,
                          reml = FALSE,
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
  
  
     # gpu coordinates
     coordsGpu = vclMatrix(coordinates, type=type)  
     # box-cox
     boxcoxGpu = vclVector(boxcox, type=type)
  
     # prepare params, make sure variance=1 in params
     params[,"variance"]=1 
     params0 = geostatsp::fillParam(params)
     paramsGpu = vclMatrix(cbind(params0, matrix(0, nrow(params0), 22-ncol(params0))),type=type)
     varMat <- vclMatrix(0, Nobs*NparamPerIter, Nobs, type=type)
     cholDiagMat <- vclMatrix(0, NparamPerIter, Nobs, type=type)
     ssqY <- vclMatrix(0, Nparam, Ndata, type=type)
     XVYXVX <- vclMatrix(0, Nparam * Ncov, ncol(yx), type=type)
     ssqBetahat <- vclMatrix(0, Nparam, Ndata, type=type)
     #ssqBeta <- vclMatrix(0, Nparam, Ndata, type=type)
     detVar <- vclVector(0, Nparam,type=type)
     detReml <- vclVector(0, Nparam, type=type)
     jacobian <- vclVector(0, Ndata, type=type)   
     ssqYX <- vclMatrix(0, ncol(yx) * NparamPerIter, ncol(yx), type=type)
     ssqYXcopy <- vclMatrix(0, ncol(yx) * NparamPerIter, ncol(yx), type=type)
     LinvYX <- vclMatrix(0, Nobs * NparamPerIter, ncol(yx), type=type)
     QinvSsqYx <- vclMatrix(0, NparamPerIter*Ncov, Ndata, type = type)
     cholXVXdiag <- vclMatrix(0, NparamPerIter, Ncov, type=type)
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
    verbose=1,  #15
    ssqYX, #
    ssqYXcopy,  #new
    LinvYX,  #18
    QinvSsqYx, 
    cholXVXdiag, #20
    varMat,        #21     Vbatch
    cholDiagMat)
  
  # resid^T V^(-1) resid, resid = Y - X betahat = ssqResidual
  ssqResidual <- ssqY - ssqBetahat
  #params0[which(is.na(as.vector(detVar))),]
  
  
  if(gpuElementsOnly==FALSE){
  if(reml== FALSE){ 
    # ml
    gpuRandom:::matrix_vector_sumBackend(Nobs*log(ssqResidual/Nobs),
                                         detVar,
                                         jacobian,  
                                         Nobs + Nobs*log(2*pi),
                                         minusTwoLogLik,
                                         Nglobal)

  }else if(reml==TRUE){
    # remlpro
    # minusTwoLogLik= (n-p)*log(ssqResidual/(n-p)) + logD + logP + jacobian + n*log(2*pi) + n-p    
    gpuRandom:::matrix_vector_sumBackend((Nobs-Ncov)*log(ssqResidual/(Nobs-Ncov)),
                                         detVar+detReml,
                                         jacobian,  
                                         Nobs*log(2*pi)+Nobs-Ncov,
                                         minusTwoLogLik,
                                         Nglobal)

  }
    
      LogLikcpu <- as.matrix(-0.5*minusTwoLogLik)
      colnames(LogLikcpu) <- paste(c('boxcox'), round(boxcox, digits = 3) ,sep = '')
  }



  if(gpuElementsOnly==FALSE){
     Output <- list(LogLik=LogLikcpu,
                    minusTwoLogLikgpu = minusTwoLogLik,
                    params = params0,
                    Nobs = Nobs,
                    Ncov = Ncov,
                    Ndata = Ndata,
                    Nparam = Nparam,
                    boxcox = boxcox,
                    jacobian = as.vector(jacobian),
                    detVar = as.vector(detVar),   
                    detReml = as.vector(detReml),   
                    ssqY = as.matrix(ssqY),    
                    XVYXVX = as.matrix(XVYXVX),
                    ssqBetahat = as.matrix(ssqBetahat),
                    ssqResidual = as.matrix(ssqResidual)
                    )
  }else{
    Output <- list(params = params0,
                   Nobs = Nobs,
                   Ncov = Ncov,
                   Ndata = Ndata,
                   Nparam = Nparam,
                   boxcox = boxcox,
                   jacobian = as.vector(jacobian),
                   detVar = as.vector(detVar),   
                   detReml = as.vector(detReml),   
                   ssqY = as.matrix(ssqY),    
                   XVYXVX = as.matrix(XVYXVX),
                   ssqBetahat = as.matrix(ssqBetahat),
                   ssqResidual = as.matrix(ssqResidual)) 
     }
     
     
     Output
     

}

