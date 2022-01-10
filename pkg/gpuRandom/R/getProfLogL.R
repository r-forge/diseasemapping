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
  
  
    if(0 %in% boxcox){
      boxcox = c(1,0, setdiff(boxcox, c(1,0)))
    }else{
      # will always be c(1,.......)
      boxcox = c(1, setdiff(boxcox, 1))
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
      # remove NAs
      selected_rows <- which(is.na(as.vector(detVar)))
      if(length(selected_rows)==0){
        paramsRenew <- params0
        detVar2 <- as.vector(detVar)
        detReml2 <- as.vector(detReml)
        ssqY2 <- as.matrix(ssqY)
        ssqBetahat2 = as.matrix(ssqBetahat)
        ssqResidual2 = as.matrix(ssqResidual)
        XVYXVX2 <- as.matrix(XVYXVX)
      }else{
      Nparam = Nparam - length(selected_rows)
      paramsRenew <- params0[-selected_rows,]
      LogLikcpu <- LogLikcpu[-selected_rows,] 
      detVar2 <- as.vector(detVar)[-selected_rows]
      detReml2 <- as.vector(detReml)[-selected_rows]
      ssqY2 <- as.matrix(ssqY)[-selected_rows,]
      ssqBetahat2 = as.matrix(ssqBetahat)[-selected_rows,]
      ssqResidual2 = as.matrix(ssqResidual)[-selected_rows,]
      
      XVYXVX2 <- as.matrix(XVYXVX)
      for (j in 1:length(selected_rows)){
        a<-c((selected_rows[j]-1)*Ncov+1, selected_rows[j]*Ncov)
        XVYXVX2 <- XVYXVX2[-a,   ]
      }
      }
  }



  if(gpuElementsOnly==FALSE){
     Output <- list(LogLik=LogLikcpu,
                    minusTwoLogLikgpu = minusTwoLogLik,
                    paramsRenew = paramsRenew,
                    Infindex = selected_rows,
                    Nobs = Nobs,
                    Ncov = Ncov,
                    Ndata = Ndata,
                    Nparam = Nparam,
                    boxcox = boxcox,
                    jacobian = as.vector(jacobian),
                    detVar = detVar2,   
                    detReml = detReml2,   
                    ssqY = ssqY2,    
                    XVYXVX = XVYXVX2,
                    ssqBetahat = ssqBetahat2,
                    ssqResidual = ssqResidual2,
                    predictors = colnames(covariates)
                    )
     
  }else{
    Output <- list(paramsRenew = paramsRenew,  # deleted Na
                   Infindex = selected_rows,
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
                   ssqResidual = as.matrix(ssqResidual),
                   predictors = colnames(covariates)) 
     }
     
     
     Output
     

  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
 

#' @export
 prof1dCov <- function(LogLik,  # cpu matrix
                       XVYXVX,  # cpu matrix
                       ssqResidual,  # cpu matrix
                       paramToEstimate = c("range","shape","nugget","anisoAngleDegrees", "anisoRatio", "boxcox"),
                       cilevel=0.95,  # decimal
                       params, # cpu matrix, 
                       boxcox,  # boxcox vallues, consistent with other functions
                       Ndata,
                       Nobs,
                       Ncov,
                       reml, 
                       predictors,  # character string
                       verbose=FALSE){
   
    
    
    
   maximum <- max(LogLik)
   breaks = maximum - qchisq(cilevel,  df = 1)/2
   
   ############## output matrix ####################
   Table <- matrix(NA, nrow=length(union(paramToEstimate, 'boxcox'))+Ncov+1, ncol=3)
   rownames(Table) <-  c(predictors, "sdSpatial", union(paramToEstimate, 'boxcox'))
   colnames(Table) <-  c("estimate", paste(c('lower', 'upper'), cilevel*100, 'ci', sep = ''))
   
   
   # myplots <- vector("list", length(paramToEstimate))
   # names(myplots) <- paramToEstimate
   if(length(paramToEstimate)==2){
   par(mfrow = c(1, 2))
   }else if(length(paramToEstimate)==3 | length(paramToEstimate)==4){
     par(mfrow = c(2, 2))
   }else if(length(paramToEstimate)==5 | length(paramToEstimate)==6){
     par(mfrow = c(2, 3))
   }else if(length(paramToEstimate)==7 | length(paramToEstimate)==8){
     par(mfrow = c(2, 2))
     par(mar = c(2.5, 3.5, 2, 0.5))
     par(mgp = c(1.5, 0.5, 0))
     par(oma = c(3, 0, 3, 0))
   }
  
   
   index <- which(LogLik == max(LogLik, na.rm = TRUE), arr.ind = TRUE)
   #################sigma hat#########################
   if(reml==FALSE)  {
     Table["sdSpatial",1] <- sqrt(ssqResidual[index[1],index[2]]/Nobs)
   }else{         
     Table["sdSpatial",1] <- sqrt(ssqResidual[index[1],index[2]]/(Nobs - Ncov))
   }
   
   params <- cbind(sqrt(params[,"nugget"]) * Table["sdSpatial",1], params)
   colnames(params)[1] <- 'sdNugget'
   
   
   ############### profile for covariance parameters #####################
   if('range' %in% paramToEstimate){
     # get profile log-lik for range values
     result = data.table::as.data.table(cbind(LogLik, params[,"range"]))
     colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 3) ,sep = ''), "range")
     profileLogLik <- result[, .(profile=max(.SD)), by=range]
     f1 <- approxfun(profileLogLik$range, profileLogLik$profile-breaks)  
     plot(profileLogLik$range, profileLogLik$profile-breaks, ylab= "proLogL-breaks", xlab='range', cex=0.5)
     # plot(profileLogLik$range, profileLogLik$profile)
     curve(f1(x), add = TRUE, col = 2, n = 1001)   #the number of x values at which to evaluate
     abline(h =0, lty = 2)
     #myplots[['range']] <- plot.range
     # rangeresults <- optim(26000, f1, method = "L-BFGS-B",lower = 20000, upper = 240000, hessian = FALSE, control=list(fnscale=-1) )
     lower = min(profileLogLik$range)
     upper = max(profileLogLik$range)
     MLE <- optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum
     # maxvalue <- rangeresults$objective
     # breaks = maxvalue - qchisq(cilevel,  df = 1)/2
     ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
     abline(v =ci[1], lty = 2)
     abline(v =ci[2], lty = 2)
     if(length(ci)==1){
       if( ci > MLE){
         ci <- c(lower, ci)
         message("did not find lower ci for range")
       }else{
         ci <- c(ci, upper)
         message("did not find upper ci for range")}
     }
     
     if(length(ci)==0 | length(ci)>2){
       warning("require a better param matrix")
       ci <- c(NA, NA)
     }
     Table["range",] <- c(MLE, ci)
     
   }
   
   
   
   if('shape' %in% paramToEstimate){
     result = as.data.table(cbind(LogLik, params[,"shape"]))
     colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 3) ,sep = ''), "shape")
     profileLogLik <- result[, .(profile=max(.SD)), by=shape]
     plot(profileLogLik$shape, profileLogLik$profile-breaks, ylab= "proLogL-breaks", xlab="shape", cex=0.5)
     f1 <- approxfun(profileLogLik$shape, profileLogLik$profile-breaks)  
     #f1 <- splinefun(profileLogLik$shape, profileLogLik$profile-breaks, method = "monoH.FC")
     curve(f1, add = TRUE, col = 2) 
     abline(h =0, lty = 2)
     #myplots[['shape']] <- plot.shape
     # shaperesults <- optim(0.1, f1, method = "L-BFGS-B",lower = 0.1, upper = 1.5, hessian = FALSE, control=list(fnscale=-1) )
     lower = min(profileLogLik$shape)
     upper = max(profileLogLik$shape)
     MLE <- optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum
     # maxvalue <- shaperesults$objective
     # breaks = maxvalue - qchisq(cilevel,  df = 1)/2
     ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
     abline(v =ci[1], lty = 2)
     abline(v =ci[2], lty = 2)
     if(length(ci)==1){
       if( ci > MLE){
         ci <- c(lower, ci)
         message("did not find lower ci for shape")
       }else{
         ci <- c(ci, upper)
         message("did not find upper ci for shape")}
     }
     
     if(length(ci)==0 | length(ci)>2){
       warning("require a better param matrix")
       ci <- c(NA, NA)
     }
     Table["shape",] <- c(MLE, ci)
     
   }       
   
   
   
   if('sdNugget' %in% paramToEstimate){
     result = as.data.table(cbind(LogLik, params[,"sdNugget"]))
     colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 2) ,sep = ''), "sdNugget")
     profileLogLik <-result[, .(profile=max(.SD)), by=sdNugget]
     plot(profileLogLik$sdNugget, profileLogLik$profile-breaks, ylab= "proLogL-breaks", xlab=bquote(tau), cex=0.5)
     f1 <- approxfun(profileLogLik$sdNugget, profileLogLik$profile-breaks)  
     # f1 <- splinefun(profileLogLik$sdNugget, profileLogLik$profile-breaks, method = "monoH.FC")
     curve(f1, add = TRUE, col = 2, n=1001) 
     abline(h =0, lty = 2)
     #myplots[['sdNugget']] <- plot.sdNugget
     # sdNuggetresults <- optim(0.1, f1, method = "L-BFGS-B",lower = 0.1, upper = 1.5, hessian = FALSE, control=list(fnscale=-1) )
     lower = min(profileLogLik$sdNugget)
     upper = max(profileLogLik$sdNugget)
     MLE <-  optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum
     # maxvalue <- sdNuggetresults$objective
     # breaks = maxvalue - qchisq(cilevel,  df = 1)/2
     ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
     abline(v =ci[1], lty = 2)
     abline(v =ci[2], lty = 2)
     if(length(ci)==1){
       if( ci > MLE){
         ci <- c(lower, ci)
         message("did not find lower ci for sdNugget")
       }else{
         ci <- c(ci, upper)
         message("did not find upper ci for sdNugget")}
     }
     
     if(length(ci)==0 | length(ci)>2){
       warning("require a better param matrix")
       ci <- c(NA, NA)
     }
     Table["sdNugget",] <- c(MLE, ci)
     
   }
   
   
   
   if('nugget' %in% paramToEstimate){
     result = as.data.table(cbind(LogLik, params[,"nugget"]))
     colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 2) ,sep = ''), "nugget")
     profileLogLik <-result[, .(profile=max(.SD)), by=nugget]
     plot(profileLogLik$nugget, profileLogLik$profile-breaks, ylab= "proLogL-breaks", xlab=expression(nu^2), cex=0.5)
     f1 <- approxfun(profileLogLik$nugget, profileLogLik$profile-breaks)  
     # f1 <- splinefun(profileLogLik$nugget, profileLogLik$profile-breaks, method = "monoH.FC")
     curve(f1, add = TRUE, col = 2, n=1001) 
     abline(h =0, lty = 2)
     #myplots[['nugget']] <- plot.nugget
     # nuggetresults <- optim(0.1, f1, method = "L-BFGS-B",lower = 0.1, upper = 1.5, hessian = FALSE, control=list(fnscale=-1) )
     lower = min(profileLogLik$nugget)
     upper = max(profileLogLik$nugget)
     MLE <-  optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum
     # maxvalue <- nuggetresults$objective
     # breaks = maxvalue - qchisq(cilevel,  df = 1)/2
     ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
     abline(v =ci[1], lty = 2)
     abline(v =ci[2], lty = 2)
     if(length(ci)==1){
       if( ci > MLE){
         ci <- c(lower, ci)
         message("did not find lower ci for nugget")
       }else{
         ci <- c(ci, upper)
         message("did not find upper ci for nugget")}
     }
     
     if(length(ci)==0 | length(ci)>2){
       warning("require a better param matrix")
       ci <- c(NA, NA)
     }
     Table["nugget",] <- c(MLE, ci)
   }
   
   
   
   
   if('anisoRatio' %in% paramToEstimate){
     result = as.data.table(cbind(LogLik, params[,"anisoRatio"]))
     colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 2) ,sep = ''), "anisoRatio")
     profileLogLik <- result[, .(profile=max(.SD)), by=anisoRatio]
     plot(profileLogLik$anisoRatio, profileLogLik$profile-breaks,ylab= "proLogL-breaks", xlab='anisoRatio', cex=0.5)
     f1 <- approxfun(profileLogLik$anisoRatio, profileLogLik$profile-breaks)  
     # f1 <- splinefun(profileLogLik$anisoRatio, profileLogLik$profile-breaks, method = "monoH.FC")
     curve(f1, add = TRUE, col = 2, n=1001) 
     abline(h =0, lty = 2)
     #myplots[['anisoRatio']] <- plot.anisoRatio
     # anisoRatioresults <- optim(0.1, f1, method = "L-BFGS-B",lower = 0.1, upper = 1.5, hessian = FALSE, control=list(fnscale=-1) )
     lower = min(profileLogLik$anisoRatio)
     upper = max(profileLogLik$anisoRatio)
     MLE <- round(optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum,digits=4)
     # maxvalue <- anisoRatioresults$objective
     # breaks = maxvalue - qchisq(cilevel,  df = 1)/2
     ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
     abline(v =ci[1], lty = 2)
     abline(v =ci[2], lty = 2)
     if(length(ci)==1){
       if( ci > MLE){
         ci <- c(lower, ci)
         message("did not find lower ci for anisoRatio")
       }else{
         ci <- c(ci, upper)
         message("did not find upper ci for anisoRatio")}
     }
     
     if(length(ci)==0 | length(ci)>2){
       warning("require a better param matrix")
       ci <- c(NA, NA)
     }
     Table["anisoRatio",] <- c(MLE, ci)
     
   }
   
   
   
   if('anisoAngleRadians' %in% paramToEstimate){
     result = as.data.table(cbind(LogLik, params[,"anisoAngleRadians"]))
     colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 2) ,sep = ''), "anisoAngleRadians")
     profileLogLik <- result[, .(profile=max(.SD)), by=anisoAngleRadians]
     plot(profileLogLik$anisoAngleRadians, profileLogLik$profile-breaks, ylab= "proLogL-breaks", xlab='anisoAngleRadians', cex=0.5)
     f1 <- approxfun(profileLogLik$anisoAngleRadians, profileLogLik$profile-breaks)
     # f1 <- splinefun(profileLogLik$anisoAngleRadians, profileLogLik$profile-breaks, method = "monoH.FC")
     curve(f1, add = TRUE, col = 2, n=1001) 
     abline(h =0, lty = 2)
     #myplots[['anisoAngleRadians']] <- plot.Degrees
     # anisoAngleRadiansresults <- optim(0.1, f1, method = "L-BFGS-B",lower = 0.1, upper = 1.5, hessian = FALSE, control=list(fnscale=-1) )
     lower = min(profileLogLik$anisoAngleRadians)
     upper = max(profileLogLik$anisoAngleRadians)
     MLE <- round(optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum,digits=4)
     # maxvalue <- anisoAngleRadiansresults$objective
     # breaks = maxvalue - qchisq(cilevel,  df = 1)/2
     ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
     abline(v =ci[1], lty = 2)
     abline(v =ci[2], lty = 2)
     if(length(ci)==1){
       if( ci > MLE){
         ci <- c(lower, ci)
         message("did not find lower ci")
       }else{
         ci <- c(ci, upper)
         message("did not find upper ci")}
     }
     
     if(length(ci)==0 | length(ci)>2){
       warning("require a better param matrix")
       ci <- c(NA, NA)
     }
     Table["anisoAngleRadians",] <- c(MLE, ci)
     
     
   }
   
   
   
   
   
   
   
   if('anisoAngleDegrees' %in% paramToEstimate){
     result = as.data.table(cbind(LogLik, params[,"anisoAngleDegrees"]))
     colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 2) ,sep = ''), "anisoAngleDegrees")
     profileLogLik <- result[, .(profile=max(.SD)), by=anisoAngleDegrees]
     plot(profileLogLik$anisoAngleDegrees, profileLogLik$profile-breaks, ylab= "proLogL-breaks", xlab='anisoAngleDegrees', cex=0.5)
     f1 <- approxfun(profileLogLik$anisoAngleDegrees, profileLogLik$profile-breaks)
     # f1 <- splinefun(profileLogLik$anisoAngleDegrees, profileLogLik$profile-breaks, method = "monoH.FC")
     curve(f1, add = TRUE, col = 2, n=1001) 
     abline(h =0, lty = 2)
     #myplots[['anisoAngleDegrees']] <- plot.Degrees
     # anisoAngleDegreesresults <- optim(0.1, f1, method = "L-BFGS-B",lower = 0.1, upper = 1.5, hessian = FALSE, control=list(fnscale=-1) )
     lower = min(profileLogLik$anisoAngleDegrees)
     upper = max(profileLogLik$anisoAngleDegrees)
     MLE <- round(optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum,digits=4)
     # maxvalue <- anisoAngleDegreesresults$objective
     # breaks = maxvalue - qchisq(cilevel,  df = 1)/2
     ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
     abline(v =ci[1], lty = 2)
     abline(v =ci[2], lty = 2)
     if(length(ci)==1){
       if( ci > MLE){
         ci <- c(lower, ci)
         message("did not find lower ci")
       }else{
         ci <- c(ci, upper)
         message("did not find upper ci")}
     }
     
     if(length(ci)==0 | length(ci)>2){
       warning("require a better param matrix")
       ci <- c(NA, NA)
     }
     Table["anisoAngleDegrees",] <- c(MLE, ci)
     
     
   }
   
   
   
   
   ###############lambda hat#####################
   if(('boxcox'%in% paramToEstimate)  & length(boxcox)>5 ){
     likForboxcox = cbind(boxcox, apply(LogLik, 2,  max) )
     f1 <- approxfun(likForboxcox[,1], likForboxcox[,2]-breaks)
     # f1 <- splinefun(likForboxcox[,1], likForboxcox[,2]-breaks, method = "monoH.FC")
     plot(likForboxcox[,1], likForboxcox[,2]-breaks, ylab= "proLogL-breaks", xlab='boxcox', cex=0.5)
     curve(f1(x), add = TRUE, col = 2, n = 1001)   #the number of x values at which to evaluate
     abline(h =0, lty = 2)
     #myplots[['boxcox']] <- plot.boxcox
     # rangeresults <- optim(26000, f1, method = "L-BFGS-B",lower = 20000, upper = 240000, hessian = FALSE, control=list(fnscale=-1) )
     lower = min(boxcox)
     upper = max(boxcox)
     MLE <- optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum
     # maxvalue <- rangeresults$objective
     # breaks = maxvalue - qchisq(cilevel,  df = 1)/2
     ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
     abline(v =ci[1], lty = 2)
     abline(v =ci[2], lty = 2)
     if(length(ci)==1){
       if( ci > MLE){
         ci <- c(lower, ci)
         message("did not find lower ci for boxcox")
       }else{
         ci <- c(ci, upper)
         message("did not find upper ci for boxcox")}
     }else if(length(ci)>2){
       warning("error in param matrix")
       ci <- c(NA, NA)
     }
     Table["boxcox",] <- c(MLE, ci)
   }else if(is.element('boxcox',paramToEstimate)  & length(boxcox) <= 5){
     message("boxcox: not enough values for interpolation!")
   }else{
     Table["boxcox",1] <- boxcox[index[2]]
   }
   
   NcolTotal = Ndata + Ncov
   
   ###############betahat#####################
   Betahat <- matrix(0, nrow=Ncov, ncol=Ndata)
   a<-c( ((index[1]-1)*Ncov+1) : (index[1]*Ncov) )
   mat <- XVYXVX[a,((Ndata+1):NcolTotal)]
   mat[upper.tri(mat)] <- mat[lower.tri(mat)]
   Betahat <- solve(mat) %*% XVYXVX[a,index[2]]
   
   Table[predictors, 1] <- Betahat

   
   
   
   Output <- list(summary = Table,
                  breaks = breaks
   )

   Output
 }      

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
#' @export
  likfitLgmCov <- function(spatialmodel,
                           params=NULL, # CPU matrix for now, users need to provide proper parameters given their specific need
                           boxcox,  # boxcox is always estimated
                           paramToEstimate, #variance and regression parameters are always estimated if not given,
                           cilevel=0.95,  # decimal
                           type = c("float", "double")[1+gpuInfo()$double_support],
                           reml=FALSE, 
                           NparamPerIter,
                           Nglobal,
                           Nlocal,
                           NlocalCache,
                           verbose=FALSE){


             data = spatialmodel$data
             formula = spatialmodel$model$formula
             coordinates = spatialmodel$data@coords
             
             if(isTRUE(params==NULL) & isTRUE(boxcox==NULL)){
               PaM <- ParamsFromLgm(spatialmodel$summary,      
                             spatialmodel$optim$mle,
                             paramToEstimate)
               ParamList <- PaM$ParamList
               params <- PaM$params
               boxcox <- PaM$boxcox
             }
             
             if(isTRUE(params==NULL) & isTRUE(boxcox!=NULL)){
               PaM <- ParamsFromLgm(spatialmodel$summary,      
                                    spatialmodel$optim$mle,
                                    paramToEstimate)
               ParamList <- PaM$ParamList
               params <- PaM$params
             }
             
             if(isTRUE(params!=NULL) & isTRUE(boxcox==NULL)){
              stop("require boxcox values to proceed")
             }
             
             result1 <- getProfLogL(data=data,
                                    formula=formula, 
                                    coordinates=coordinates,
                                    params=params,  # CPU matrix 
                                    boxcox=boxcox,  # boxcox is always estimated
                                    type = type,
                                    NparamPerIter =NparamPerIter,
                                    gpuElementsOnly = FALSE,
                                    reml = reml,
                                    Nglobal, Nlocal, NlocalCache, verbose=FALSE)
             

             
             result2 <- prof1dCov(LogLik = result1$LogLik,  # cpu matrix
                                  XVYXVX = result1$XVYXVX,  # cpu matrix
                                  ssqResidual = result1$ssqResidual,  # cpu matrix
                                  paramToEstimate = paramToEstimate,
                                  cilevel=cilevel,  # decimal
                                  params = result1$paramsRenew, # cpu matrix, 
                                  boxcox = result1$boxcox,  # boxcox vallues, consistent with other functions
                                  Ndata = result1$Ndata,
                                  Nobs = result1$Nobs,
                                  Ncov = result1$Ncov,
                                  reml = reml, 
                                  predictors = result1$predictors,
                                  verbose=FALSE)
             
              

             if(isTRUE(params==NULL)){
               Output <- list(summary = result2$summary,
                              breaks = result2$breaks,
                              LogLik = result1$LogLik,
                              minusTwoLogLikgpu = result1$minusTwoLogLikgpu,
                              params = result1$params,
                              ParamList = ParamList,
                              boxcox = result1$boxcox,
                              Nobs = result1$Nobs,
                              Ncov = result1$Ncov,
                              Ndata = result1$Ndata,
                              Nparam = result1$Nparam,
                              jacobian = result1$jacobian,
                              detVar = result1$detVar,   
                              detReml = result1$detReml,   
                              ssqY = result1$ssqY,    
                              XVYXVX = result1$XVYXVX,
                              ssqBetahat = result1$ssqBetahat,
                              ssqResidual = result1$ssqResidual)   
               
             }else{
             Output <- list(summary = result2$summary,
                            breaks = result2$breaks,
                            LogLik = result1$LogLik,
                            minusTwoLogLikgpu = result1$minusTwoLogLikgpu,
                            boxcox = result1$boxcox,
                            Nobs = result1$Nobs,
                            Ncov = result1$Ncov,
                            Ndata = result1$Ndata,
                            Nparam = result1$Nparam,
                            jacobian = result1$jacobian,
                            detVar = result1$detVar,   
                            detReml = result1$detReml,   
                            ssqY = result1$ssqY,    
                            XVYXVX = result1$XVYXVX,
                            ssqBetahat = result1$ssqBetahat,
                            ssqResidual = result1$ssqResidual)
             }     
             
             
             Output
             
  }           
             
             
             

             
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 