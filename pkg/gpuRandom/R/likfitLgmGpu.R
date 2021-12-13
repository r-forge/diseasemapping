#' @title Estimate profile Log-likelihood for covariance parameters and lambda
#' @import data.table
#' @importFrom rootSolve uniroot.all
#' @useDynLib gpuRandom
#' @export



# betahat and sigmahat 
# profile log-likelihood for each covariance parameters + lambda
likfitLgmGpu <- function(data,
                         formula, 
                         coordinates,
                         params, # CPU matrix for now, users need to provide proper parameters given their specific need
                         #fixVariance = FALSE, 
                         boxcox,  # boxcox is always estimated
                         paramToEstimate = c("range","shape","nugget","anisoAngleDegrees", "anisoRatio", "boxcox"), #variance and regression parameters are always estimated if not given,
                         cilevel=0.95,  # decimal
                         type = c("float", "double")[1+gpuInfo()$double_support],
                         reml=FALSE, 
                         minustwotimes=FALSE,
                         NparamPerIter,
                         Nglobal,
                         Nlocal,
                         NlocalCache,
                         verbose=FALSE){
  
  # will always be c(1,0,.......)
  boxcox = c(1, 0, setdiff(boxcox, c(0,1)))
  
  if(length(boxcox)<=3 & is.element('boxcox',paramToEstimate))
  {
    warning("need more params to estimate boxcox")
  }
  
  
  
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
  params0 = geostatsp::fillParam(params)
  params0[,"variance"]=1    
  paramsGpu = vclMatrix(cbind(params0, matrix(0, nrow(params0), 22-ncol(params0))),type=type)
  
  varMat = vclMatrix(0, Nobs*NparamPerIter, Nobs, type=type)
  cholDiagMat = vclMatrix(0, NparamPerIter, Nobs, type=type)
  ssqY <- vclMatrix(0, Nparam, Ndata, type=type)
  XVYXVX = vclMatrix(0, Nparam * Ncov, ncol(yx), type=type)
  ssqBetahat <- vclMatrix(0, Nparam, Ndata, type=type)
  ssqBeta <- vclMatrix(0, Nparam, Ndata, type=type)
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
  # any(is.na(as.matrix(log(ssqResidual/Nobs))))
   # as.matrix(cholDiagMat)[]
   any(is.nan(as.matrix(log(ssqResidual/Nobs))))
   any(is.nan(as.vector(detVar)))
   any(is.na(as.matrix(varMat)))
   any(is.nan(as.matrix(varMat)))
   any(is.na(as.matrix(ssqYX)))
   any(is.na(as.matrix(ssqY)))
   any(is.na(as.matrix(ssqYXcopy)))
   any(is.na(as.matrix(ssqBetahat)))
   any(is.na(as.vector(detReml)))
   any(is.na(as.matrix(cholXVXdiag)))
   any(is.na(as.matrix(ssqResidual)))
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
  # 
  # params0[which(is.na(as.vector(detVar))),]
  # 
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
  breaks = maximum - qchisq(cilevel,  df = 1)/2
  
  ############## output matrix ####################
  Table <- matrix(NA, nrow=length(union(paramToEstimate, 'boxcox'))+Ncov+1, ncol=3)
  rownames(Table) <-  c("intercept", paste(c('betahat'), seq_len(Ncov-1),sep = ''), "sdSpatial", union(paramToEstimate, 'boxcox'))
  colnames(Table) <-  c("estimate", paste(c('lower', 'upper'), cilevel*100, 'ci', sep = ''))
  
  

  
  ############### profile for covariance parameters #####################
  if(is.element('range',paramToEstimate)){
    # get profile log-lik for range values
    result = data.table::as.data.table(cbind(LogLikcpu, params0[,"range"]))
    colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 2) ,sep = ''), "range")
    profileLogLik <- result[, .(profile=max(.SD)), by=range]
    f1 <- approxfun(profileLogLik$range, profileLogLik$profile-breaks)  
    # f1 <- splinefun(profileLogLik$range, profileLogLik$profile-breaks, method = "monoH.FC")
    # plot(profileLogLik$range, profileLogLik$profile-breaks)
    # plot(profileLogLik$range, profileLogLik$profile)
    # curve(f1(x), add = TRUE, col = 2, n = 1001)   #the number of x values at which to evaluate
    # abline(h =0, lty = 2)
    # rangeresults <- optim(26000, f1, method = "L-BFGS-B",lower = 20000, upper = 240000, hessian = FALSE, control=list(fnscale=-1) )
    lower = min(profileLogLik$range)
    upper = max(profileLogLik$range)
    MLE <- optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum
    # maxvalue <- rangeresults$objective
    # breaks = maxvalue - qchisq(cilevel,  df = 1)/2
    ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
    if(length(ci)==1){
      if( ci > MLE){
        ci <- c(lower, ci)
        warning("did not find lower ci for range, require more params")
      }else{
        ci <- c(ci, upper)
        warning("did not find upper ci for range, require more params")}
    }
    
    if(length(ci)==0 | length(ci)>2){
      warning("require a better param matrix")
      ci <- c(NA, NA)
    }
    Table["range",] <- c(MLE, ci)
    # MLE <- profileLogLik$range[which.max(profileLogLik$profile)]
    # plot(profileLogLik$range, profileLogLik$profile)
    # f1<-approxfun(profileLogLik$range, profileLogLik$profile-breaks)  
    # optimization <- optimize(f1, interval = c(min(profileLogLik$range), max(profileLogLik$range)),maximum = TRUE)
    # MLE <- optimization$maximum 
    # uniroot(f1, lower = 20000, upper = 80000)$root
    # abline(v =ci[1], lty = 2)
    # abline(v =ci[2], lty = 2)
    # points(All, y = rep(0, length(All)), pch = 15, cex =0.7 )
    # All2 <-c(afLeft(breaks95), afRight(breaks95))
    # abline(v =All2[1], lty = 2)
    # abline(v =All2[2], lty = 2)
    # points(All2, y = rep(0, length(All2)), col = "red", cex = 0.7)
    
    
    # leftOfMax = profileLogLik$range < MLE
    # if(length(which(leftOfMax)) <=2 | length(which(!leftOfMax)) <=2){
    #   ci95 = c(NA,NA)
    #   print("Not enough data for CI calculation")
    # }else{
    #   afLeft <- approxfun(profileLogLik$profile[leftOfMax], profileLogLik$range[leftOfMax])   
    #   afRight <- approxfun(profileLogLik$profile[!leftOfMax], profileLogLik$range[!leftOfMax])   
    #   # plot(profileLogLik$profile[leftOfMax], profileLogLik$range[leftOfMax])
    #   # plot(profileLogLik$profile[!leftOfMax], profileLogLik$range[!leftOfMax])
    #   # curve(afLeft, add=TRUE)
    #   # curve(afRight, add=TRUE)
    #   
    #  
    #   ci95= c(afLeft(breaks95), afRight(breaks95))
    #   if (any(is.na(ci95))){
    #     warning("'range' scope too small for 95% ci")
    #   }
    # }

  }
  
  
  
  
  
  if(is.element('shape',paramToEstimate)){
    result = as.data.table(cbind(LogLikcpu, params0[,"shape"]))
    colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 2) ,sep = ''), "shape")
    profileLogLik <- result[, .(profile=max(.SD)), by=shape]
    # plot(profileLogLik$shape, profileLogLik$profile-breaks)
    f1 <- approxfun(profileLogLik$shape, profileLogLik$profile-breaks)  
    #f1 <- splinefun(profileLogLik$shape, profileLogLik$profile-breaks, method = "monoH.FC")
    # curve(f1, add = TRUE, col = 2) 
    # abline(h =0, lty = 2)
    # shaperesults <- optim(0.1, f1, method = "L-BFGS-B",lower = 0.1, upper = 1.5, hessian = FALSE, control=list(fnscale=-1) )
    lower = min(profileLogLik$shape)
    upper = max(profileLogLik$shape)
    MLE <- optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum
    # maxvalue <- shaperesults$objective
    # breaks = maxvalue - qchisq(cilevel,  df = 1)/2
    ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
    if(length(ci)==1){
      if( ci > MLE){
        ci <- c(lower, ci)
        warning("did not find lower ci for shape, require more params")
      }else{
        ci <- c(ci, upper)
        warning("did not find upper ci for shape, require more params")}
    }
    
    if(length(ci)==0 | length(ci)>2){
      warning("require a better param matrix")
      ci <- c(NA, NA)
    }
    Table["shape",] <- c(MLE, ci)

    # result = as.data.table(cbind(LogLikcpu, params0[,"shape"]))
    # colnames(result) <- c(paste(c('boxcox'),boxcox,sep = ''), "shape")
    # profileLogLik <- result[, .(profile=max(.SD)), by=shape]
    # MLE <- profileLogLik$shape[which.max(profileLogLik$profile)]
    # #plot(profileLogLik$shape, profileLogLik$profile)
    # leftOfMax = profileLogLik$shape < profileLogLik$shape[which.max(profileLogLik$profile)]
    # if(length(which(leftOfMax)) <=2 | length(which(!leftOfMax)) <=2){
    #   ci95 = c(NA,NA)
    #   print("Not enough data for CI calculation")
    # }else{
    #   afLeft <- approxfun(profileLogLik$profile[leftOfMax], profileLogLik$shape[leftOfMax])  
    #   afRight <- approxfun(profileLogLik$profile[!leftOfMax], profileLogLik$shape[!leftOfMax])   
    #   breaks95 = maximum - qchisq(0.95,  df = 1)/2
    #   ci95= c(afLeft(breaks95), afRight(breaks95))
    # }
    # 
    # 
    # if (any(is.na(ci95)))
    #   warning("'shape' scope too small for 95% ci")
}       
  
  
  
  
  
  if(is.element('nugget',paramToEstimate)){
    result = as.data.table(cbind(LogLikcpu, params0[,"nugget"]))
    colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 2) ,sep = ''), "nugget")
    profileLogLik <-result[, .(profile=max(.SD)), by=nugget]
    # plot(profileLogLik$nugget, profileLogLik$profile-breaks)
    f1 <- approxfun(profileLogLik$nugget, profileLogLik$profile-breaks)  
    # f1 <- splinefun(profileLogLik$nugget, profileLogLik$profile-breaks, method = "monoH.FC")
    # curve(f1, add = TRUE, col = 2, n=1001) 
    # abline(h =0, lty = 2)
    # nuggetresults <- optim(0.1, f1, method = "L-BFGS-B",lower = 0.1, upper = 1.5, hessian = FALSE, control=list(fnscale=-1) )
    lower = min(profileLogLik$nugget)
    upper = max(profileLogLik$nugget)
    MLE <-  optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum
    # maxvalue <- nuggetresults$objective
    # breaks = maxvalue - qchisq(cilevel,  df = 1)/2
    ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
    if(length(ci)==1){
      if( ci > MLE){
        ci <- c(lower, ci)
        warning("did not find lower ci for nugget, require more params")
      }else{
        ci <- c(ci, upper)
        warning("did not find upper ci for nugget, require more params")}
    }
    
    if(length(ci)==0 | length(ci)>2){
      warning("require a better param matrix")
      ci <- c(NA, NA)
    }
    Table["nugget",] <- c(MLE, ci)
    
    # MLE <- profileLogLik$nugget[which.max(profileLogLik$profile)]
    # leftOfMax = profileLogLik$nugget < MLE
    # if(length(which(leftOfMax)) <=2 | length(which(!leftOfMax)) <=2){
    #   ci95 = c(NA,NA)
    #   print("Not enough data for CI calculation")
    # }else{
    #   afLeft <- approxfun(profileLogLik$profile[leftOfMax], profileLogLik$nugget[leftOfMax], rule=1)   # plot(profileLogLik$range, profileLogLik$profile)# curve(af, add=TRUE)
    #   afRight <- approxfun(profileLogLik$profile[!leftOfMax], profileLogLik$nugget[!leftOfMax])   
    #   # plot(profileLogLik$nugget, profileLogLik$profile)
    #   # curve(afRight, add=TRUE)
    #   ci95= c(afLeft(breaks95), afRight(breaks95))
    # }
    # if (any(is.na(ci95)))
    #   warning("'nugget' scope too small for 95% ci")
    
  }
  
  
  
  
  if(is.element('anisoRatio',paramToEstimate)){
    result = as.data.table(cbind(LogLikcpu, params0[,"anisoRatio"]))
    colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 2) ,sep = ''), "anisoRatio")
    profileLogLik <- result[, .(profile=max(.SD)), by=anisoRatio]
    # plot(profileLogLik$anisoRatio, profileLogLik$profile-breaks)
    f1 <- approxfun(profileLogLik$anisoRatio, profileLogLik$profile-breaks)  
    # f1 <- splinefun(profileLogLik$anisoRatio, profileLogLik$profile-breaks, method = "monoH.FC")
    # curve(f1, add = TRUE, col = 2, n=1001) 
    # abline(h =0, lty = 2)
    # anisoRatioresults <- optim(0.1, f1, method = "L-BFGS-B",lower = 0.1, upper = 1.5, hessian = FALSE, control=list(fnscale=-1) )
    lower = min(profileLogLik$anisoRatio)
    upper = max(profileLogLik$anisoRatio)
    MLE <- round(optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum,digits=4)
    # maxvalue <- anisoRatioresults$objective
    # breaks = maxvalue - qchisq(cilevel,  df = 1)/2
    ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
    if(length(ci)==1){
      if( ci > MLE){
        ci <- c(lower, ci)
        warning("did not find lower ci for anisoRatio, require more params")
      }else{
        ci <- c(ci, upper)
        warning("did not find upper ci for anisoRatio, require more params")}
    }
    
    if(length(ci)==0 | length(ci)>2){
      warning("require a better param matrix")
      ci <- c(NA, NA)
    }
    Table["anisoRatio",] <- c(MLE, ci)

    # maximum <- max(profileLogLik$profile)
    # MLE <- profileLogLik$anisoRatio[which.max(profileLogLik$profile)]
    # 
    # leftOfMax = profileLogLik$anisoRatio < MLE
    # if(length(which(leftOfMax)) <=2 | length(which(!leftOfMax)) <=2){
    #   ci95 = c(NA,NA)
    #   print("Not enough data for CI calculation")
    # }else{
    #   afLeft <- approxfun(profileLogLik$profile[leftOfMax], profileLogLik$anisoRatio[leftOfMax])   
    #   # plot(profileLogLik$anisoRatio, profileLogLik$profile)# curve(af, add=TRUE)
    #   afRight <- approxfun(profileLogLik$profile[!leftOfMax], profileLogLik$anisoRatio[!leftOfMax])   
    #   
    #   ci95= c(afLeft(breaks95), afRight(breaks95))}
    # if (any(is.na(ci95)))
    #   warning("'anisoRatio' scope too small for 95% ci")

  }
  
  
  
  if(is.element('anisoAngleDegrees',paramToEstimate)){
    result = as.data.table(cbind(LogLikcpu, params0[,"anisoAngleDegrees"]))
    colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 2) ,sep = ''), "anisoAngleDegrees")
    profileLogLik <- result[, .(profile=max(.SD)), by=anisoAngleDegrees]
    #plot(profileLogLik$anisoAngleDegrees, profileLogLik$profile-breaks)
    f1 <- approxfun(profileLogLik$anisoAngleDegrees, profileLogLik$profile-breaks)
    # f1 <- splinefun(profileLogLik$anisoAngleDegrees, profileLogLik$profile-breaks, method = "monoH.FC")
    # curve(f1, add = TRUE, col = 2, n=1001) 
    # abline(h =0, lty = 2)
    # anisoAngleDegreesresults <- optim(0.1, f1, method = "L-BFGS-B",lower = 0.1, upper = 1.5, hessian = FALSE, control=list(fnscale=-1) )
    lower = min(profileLogLik$anisoAngleDegrees)
    upper = max(profileLogLik$anisoAngleDegrees)
    MLE <- round(optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum,digits=4)
    # maxvalue <- anisoAngleDegreesresults$objective
    # breaks = maxvalue - qchisq(cilevel,  df = 1)/2
    ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
    if(length(ci)==1){
      if( ci > MLE){
        ci <- c(lower, ci)
        warning("did not find lower ci, require more params")
      }else{
        ci <- c(ci, upper)
        warning("did not find upper ci, require more params")}
    }
    
    if(length(ci)==0 | length(ci)>2){
      warning("require a better param matrix")
      ci <- c(NA, NA)
    }
    Table["anisoAngleDegrees",] <- c(MLE, ci)

    
    # MLE <- profileLogLik$anisoAngleDegrees[which.max(profileLogLik$profile)]
    # leftOfMax = profileLogLik$anisoAngleDegrees < MLE
    # if(length(which(leftOfMax)) <=2 | length(which(!leftOfMax)) <=2){
    #   ci95 = c(NA,NA)
    #   print("Not enough data for CI calculation")
    # }else{
    #   afLeft <- approxfun(profileLogLik$profile[leftOfMax], profileLogLik$anisoAngleDegrees[leftOfMax])   
    #   # plot(profileLogLik$anisoAngleDegrees, profileLogLik$profile)# curve(af, add=TRUE)
    #   afRight <- approxfun(profileLogLik$profile[!leftOfMax], profileLogLik$anisoAngleDegrees[!leftOfMax])   
    #   ci95= c(afLeft(breaks95), afRight(breaks95))}
    # if (any(is.na(ci95)))
    #   warning("'anisoAngleDegrees' scope too small for 95% ci")

  }
  
  
  
  index <- which(LogLikcpu == max(LogLikcpu, na.rm = TRUE), arr.ind = TRUE)
  
  ###############lambda hat#####################
  if(is.element('boxcox',paramToEstimate)  & length(boxcox)>3 ){
    likForboxcox = cbind(boxcox, apply(LogLikcpu, 2, max) )
    f1 <- approxfun(likForboxcox[,1], likForboxcox[,2]-breaks)
    # f1 <- splinefun(likForboxcox[,1], likForboxcox[,2]-breaks, method = "monoH.FC")
    # plot(likForboxcox[,1], likForboxcox[,2]-breaks)
    # curve(f1(x), add = TRUE, col = 2, n = 1001)   #the number of x values at which to evaluate
    # abline(h =0, lty = 2)
    # rangeresults <- optim(26000, f1, method = "L-BFGS-B",lower = 20000, upper = 240000, hessian = FALSE, control=list(fnscale=-1) )
    lower = min(boxcox)
    upper = max(boxcox)
    MLE <- optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum
    # maxvalue <- rangeresults$objective
    # breaks = maxvalue - qchisq(cilevel,  df = 1)/2
    ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
    if(length(ci)==1){
      if( ci > MLE){
        ci <- c(lower, ci)
        warning("did not find lower ci for boxcox, require more params")
      }else{
        ci <- c(ci, upper)
        warning("did not find upper ci for boxcox, require more params")}
    }else if(length(ci)>2){
      warning("error in param matrix")
      ci <- c(NA, NA)
    }
    Table["boxcox",] <- c(MLE, ci)
  }else{
    Table["boxcox",1] <- boxcox[index[2]]
  }
  


  ###############betahat#####################
  Betahat <- matrix(0, nrow=Ncov, ncol=Ndata)
  a<-c( ((index[1]-1)*Ncov+1) : (index[1]*Ncov) )
  mat <- XVYXVX[a,((Ndata+1):ncol(yx))]
  mat[upper.tri(mat)] <- mat[lower.tri(mat)]
  Betahat <- solve(mat) %*% XVYXVX[a,index[2]]
  
  Table[c("intercept", paste(c('betahat'), seq_len(Ncov-1),sep = '')),1] <- Betahat
  
  
  
  
  #################sigma hat#########################
  if(reml==FALSE)  {
    Table["sdSpatial",1] <- sqrt(ssqResidual[index[1],index[2]]/Nobs)
  }
  else if(reml==TRUE) {         
    Table["sdSpatial",1] <- sqrt(ssqResidual[index[1],index[2]]/(Nobs - Ncov))
  }
  
  
  
  Output <- list(LogLik=LogLikcpu,
                 minusTwoLogLik = minusTwoLogLik,
                 estimates = Table,
                 Nobs = Nobs,
                 Ncov = Ncov,
                 Ndata = Ndata,
                 Nparam = Nparam,
                 breaks = breaks,
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












