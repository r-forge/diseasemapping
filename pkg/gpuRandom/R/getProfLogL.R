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
     params0 = geostatsp::fillParam(params)
     params0[,"variance"]=1 
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
  
  
  
  
  
  
  
  
  

  
  
 
#' @import mgcv
#' @import data.table
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
   # if(length(paramToEstimate)==2){
   par(mfrow = c(1, 1))
   # }else if(length(paramToEstimate)==3 | length(paramToEstimate)==4){
   #   par(mfrow = c(2, 2))
   # }else if(length(paramToEstimate)==5 | length(paramToEstimate)==6){
   #   par(mfrow = c(2, 3))
   # }else if(length(paramToEstimate)==7 | length(paramToEstimate)==8){
   #   par(mfrow = c(2, 2))
   #   par(mar = c(2.5, 3.5, 1, 0.5))
   #   par(mgp = c(1.5, 0.5, 0))
   #   par(oma = c(0, 0, 3, 0))
   # }
  
   
   index <- which(LogLik == max(LogLik, na.rm = TRUE), arr.ind = TRUE)
   #################sigma hat#########################
   if(reml==FALSE)  {
     Table["sdSpatial",1] <- sqrt(ssqResidual[index[1],index[2]]/Nobs)
   }else{         
     Table["sdSpatial",1] <- sqrt(ssqResidual[index[1],index[2]]/(Nobs - Ncov))
   }
   
   
############### profile for covariance parameters #####################
######################range ########
   if('range' %in% paramToEstimate){
     result = data.table::as.data.table(cbind(LogLik, params[,"range"]))
     colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 3) ,sep = ''), "range")
     profileLogLik <- result[, .(profile=max(.SD)), by=range]
     profileLogLik[,'profile'] <- profileLogLik[,'profile'] - breaks

     datC1= geometry::convhulln(profileLogLik)
     allPoints1 = unique(as.vector(datC1))
     toTest = profileLogLik[allPoints1,]
     toTest[,'profile'] = toTest[,'profile'] - 0.01
     inHull1 = geometry::inhulln(datC1, as.matrix(toTest))
     toUse = profileLogLik[allPoints1,][inHull1,]
     
     interp1 = mgcv::gam(profile ~ s(range, k=nrow(toUse), fx=TRUE), data=toUse)
     #prof1 = data.frame(range=seq(min(toUse[,1])-0.1, max(toUse[,1])+0.1, len=1001))
     #prof1$z = predict(interp1, prof1)
     
     plot(profileLogLik$range, profileLogLik$profile, cex=.2, xlab="range", ylab="profileLogL")
     points(toTest, col='red', cex=0.6)
     points(toUse, col='blue', cex=0.6, pch=3)
     f1 <- approxfun(toUse$range, toUse$profile)
     curve(f1(x), add = TRUE, col = 'green', n = 1001)   #the number of x values at which to evaluate
     abline(h =0, lty = 2, col='red')
     lower = min(profileLogLik$range)
     upper = max(profileLogLik$range)
     MLE <- optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum
     
     ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
     abline(v =ci[1], lty = 2, col='red')
     abline(v =ci[2], lty = 2, col='red')
     if(length(ci)==1){
        if(ci > MLE){
           ci <- c(lower, ci)
           message("did not find lower ci for range")
        }else{
           ci <- c(ci, upper)
           message("did not find upper ci for range")}
     }
     
     if(length(ci)==0 | length(ci)>2){
        warning("error in params")
        ci <- c(NA, NA)
     }
     Table["range",] <- c(MLE, ci)
     ####################### log plot ##############################
     profileLogLik$logrange <- log(profileLogLik$range)
     datC1= geometry::convhulln(profileLogLik[,c('logrange','profile')])
     allPoints1 = unique(as.vector(datC1))
     toTest = profileLogLik[allPoints1,c('logrange','profile')]
     toTest[,'profile'] = toTest[,'profile'] - 0.01
     inHull1 = geometry::inhulln(datC1, as.matrix(toTest))
     toUse = profileLogLik[allPoints1,c('logrange','profile')][inHull1,]
     
     interp1 = mgcv::gam(profile ~ s(logrange, k=nrow(toUse), fx=TRUE), data=toUse)
     prof1 = data.frame(logrange=seq(min(toUse[,1])-0.1, max(toUse[,1])+0.1, len=501))
     prof1$z = predict(interp1, prof1)
     
     plot(profileLogLik$logrange, profileLogLik$profile, cex=.2, xlab="log(range)", ylab="profileLogL")
     points(toTest, col='red', cex=0.6)
     points(toUse, col='blue', cex=0.6, pch=3)
     lines(prof1$logrange, prof1$z, col='green')
     abline(h =0, lty = 2, col='red')
}
   
   
   
   
################shape ##############   
   if('shape' %in% paramToEstimate){
     result = as.data.table(cbind(LogLik, params[,"shape"]))
     colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 3) ,sep = ''), "shape")
     profileLogLik <- result[, .(profile=max(.SD)), by=shape]
     profileLogLik[,'profile'] <- profileLogLik[,'profile'] - breaks
     
     datC1= geometry::convhulln(profileLogLik)
     allPoints1 = unique(as.vector(datC1))
     toTest = profileLogLik[allPoints1,]
     toTest[,'profile'] = toTest[,'profile'] - 0.01
     inHull1 = geometry::inhulln(datC1, as.matrix(toTest))
     toUse = profileLogLik[allPoints1,][inHull1,]
     
     interp1 = mgcv::gam(profile ~ s(shape, k=nrow(toUse), fx=TRUE), data=toUse)
     prof1 = data.frame(shape=seq(min(toUse[,1])-0.1, max(toUse[,1])+0.1, len=501))
     prof1$z = predict(interp1, prof1)
     
     plot(profileLogLik$shape, profileLogLik$profile, cex=.2, xlab="shape", ylab="profileLogL", xlim=c(0,50))
     points(toTest, col='red', cex=0.6)
     points(toUse, col='blue', cex=0.6, pch=3)
     lines(prof1$shape, prof1$z, col='green')
     f1 <- approxfun(prof1$shape, prof1$z)
     #f1 <- approxfun(toUse$shape, toUse$profile)
     #curve(f1(x), add = TRUE, col = 'green', n = 1001)   #the number of x values at which to evaluate
     abline(h =0, lty = 2, col='red')
     lower = min(profileLogLik$shape)
     upper = 50
     MLE <- optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum
     
     ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
     abline(v =ci[1], lty = 2, col='red')
     abline(v =ci[2], lty = 2, col='red')
     if(length(ci)==1){
        if(ci > MLE){
           ci <- c(lower, ci)
           message("did not find lower ci for shape")
        }else{
           ci <- c(ci, upper)
           message("did not find upper ci for shape")}
     }
     
     if(length(ci)==0 | length(ci)>2){
        warning("error in params")
        ci <- c(NA, NA)
     }
     Table["shape",] <- c(MLE, ci)
     ####################### log plot ##############################
     profileLogLik$logshape <- log(profileLogLik$shape)
     datC1= geometry::convhulln(profileLogLik[,c('logshape','profile')])
     allPoints1 = unique(as.vector(datC1))
     toTest = profileLogLik[allPoints1,c('logshape','profile')]
     toTest[,'profile'] = toTest[,'profile'] - 0.01
     inHull1 = geometry::inhulln(datC1, as.matrix(toTest))
     toUse = profileLogLik[allPoints1,c('logshape','profile')][inHull1,]
     
     interp1 = mgcv::gam(profile ~ s(logshape, k=nrow(toUse), fx=TRUE), data=toUse)
     prof1 = data.frame(logshape=seq(min(toUse[,1])-0.1, max(toUse[,1])+0.1, len=501))
     prof1$z = predict(interp1, prof1)
     
     plot(profileLogLik$logshape, profileLogLik$profile, cex=.2, xlab="log(shape)", ylab="profileLogL")
     points(toTest, col='red', cex=0.6)
     points(toUse, col='blue', cex=0.6, pch=3)
     lines(prof1$logshape, prof1$z, col='green')
     abline(h =0, lty = 2, col='red')
     
   }       
   
   
   
################sd nugget ##############     
   params <- cbind(params, sqrt(params[,"nugget"]) * Table["sdSpatial",1])
   colnames(params)[ncol(params)] <- 'sdNugget'
   
   if('sdNugget' %in% paramToEstimate){
      result = data.table::as.data.table(cbind(LogLik, params[,"sdNugget"]))
      colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 3) ,sep = ''), "sdNugget")
      profileLogLik <- result[, .(profile=max(.SD)), by=sdNugget]
      profileLogLik[,'profile'] <- profileLogLik[,'profile'] - breaks
      
      datC1= geometry::convhulln(profileLogLik)
      allPoints1 = unique(as.vector(datC1))
      toTest = profileLogLik[allPoints1,]
      toTest[,'profile'] = toTest[,'profile'] - 0.01
      inHull1 = geometry::inhulln(datC1, as.matrix(toTest))
      toUse = profileLogLik[allPoints1,][inHull1,]
      toUse = toUse[order(toUse[,1])]
      interp1 = mgcv::gam(profile ~ s(sdNugget, k=nrow(toUse), fx=TRUE), data=toUse)
      prof1 = data.frame(sdNugget=seq(min(toUse[,1])-0.1, max(toUse[,1])+0.1, len=1001))
      prof1$z = predict(interp1, prof1)
      
      plot(profileLogLik$sdNugget, profileLogLik$profile, cex=.2, xlab="sdNugget", ylab="profileLogL")
      points(toTest, col='red', cex=0.6)
      points(toUse, col='blue', cex=0.6, pch=3)
      lines(prof1$sdNugget, prof1$z, col='green')
      f1 <- approxfun(prof1$sdNugget, prof1$z)
      #f1 <- approxfun(toUse$sdNugget, toUse$profile)
      #curve(f1(x), add = TRUE, col = 'green', n = 1001)   #the number of x values at which to evaluate
      abline(h =0, lty = 2, col='red')
      lower = min(profileLogLik$sdNugget)
      upper = max(profileLogLik$sdNugget)
      MLE <- optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum
      
      ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
      abline(v =ci[1], lty = 2, col='red')
      abline(v =ci[2], lty = 2, col='red')
      if(length(ci)==1){
         if(ci > MLE){
            ci <- c(lower, ci)
            message("did not find lower ci for sdNugget")
         }else{
            ci <- c(ci, upper)
            message("did not find upper ci for sdNugget")}
      }
      
      if(length(ci)==0 | length(ci)>2){
         warning("error in params")
         ci <- c(NA, NA)
      }
      Table["sdNugget",] <- c(MLE, ci)
      ####################### log plot ##############################
      # plot(log(profileLogLik$sdNugget), profileLogLik$profile, cex=.2, xlab="log(sdNugget)", ylab="profileLogL")
      # toTestLog = cbind(log(toTest[,'sdNugget']), toTest[,'profile'])
      # toUseLog = cbind(log(toUse[,'sdNugget']), toUse[,'profile'])
      # points(toTestLog, col='red', cex=0.6)
      # points(toUseLog, col='blue', cex=0.6, pch=3)
      # toUseLog <- toUseLog[order(toUseLog[,1],decreasing=FALSE)]
      # lines(toUseLog$sdNugget, toUseLog$profile, col='green')   #the number of x values at which to evaluate
      # abline(h =0, lty = 2, col='red')

   }
   
   
   
   if('nugget' %in% paramToEstimate){
      result = as.data.table(cbind(LogLik, params[,"nugget"]))
      colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 3) ,sep = ''), "nugget")
      profileLogLik <- result[, .(profile=max(.SD)), by=nugget]
      profileLogLik[,'profile'] <- profileLogLik[,'profile'] - breaks
      
      datC1= geometry::convhulln(profileLogLik)
      allPoints1 = unique(as.vector(datC1))
      toTest = profileLogLik[allPoints1,]
      toTest[,'profile'] = toTest[,'profile'] - 0.01
      inHull1 = geometry::inhulln(datC1, as.matrix(toTest))
      toUse = profileLogLik[allPoints1,][inHull1,]
      
      interp1 = mgcv::gam(profile ~ s(nugget, k=nrow(toUse), fx=TRUE), data=toUse)
      prof1 = data.frame(nugget=seq(min(toUse[,1])-0.1, max(toUse[,1])+0.1, len=501))
      prof1$z = predict(interp1, prof1)
      
      plot(profileLogLik$nugget, profileLogLik$profile, cex=.2, xlab="nugget", ylab="profileLogL")
      points(toTest, col='red', cex=0.6)
      points(toUse, col='blue', cex=0.6, pch=3)
      lines(prof1$nugget, prof1$z, col='green')
      prof1 <- prof1[prof1$nugget>0,]
      f1 <- approxfun(prof1$nugget, prof1$z)
      #f1 <- approxfun(toUse$nugget, toUse$profile)
      #curve(f1(x), add = TRUE, col = 'green', n = 1001)   #the number of x values at which to evaluate
      abline(h =0, lty = 2, col='red')
      lower = min(profileLogLik$nugget)
      upper = max(profileLogLik$nugget)
      MLE <- optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum
      
      ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
      abline(v =ci[1], lty = 2, col='red')
      abline(v =ci[2], lty = 2, col='red')
      if(length(ci)==1){
         if(ci > MLE){
            ci <- c(lower, ci)
            message("did not find lower ci for nugget")
         }else{
            ci <- c(ci, upper)
            message("did not find upper ci for nugget")}
      }
      
      if(length(ci)==0 | length(ci)>2){
         warning("error in params")
         ci <- c(NA, NA)
      }
      Table["nugget",] <- c(MLE, ci)
      ####################### log plot ##############################
      profileLogLik$lognugget <- log(profileLogLik$nugget)
      datC1= geometry::convhulln(profileLogLik[,c('lognugget','profile')])
      allPoints1 = unique(as.vector(datC1))
      toTest = profileLogLik[allPoints1,c('lognugget','profile')]
      toTest[,'profile'] = toTest[,'profile'] - 0.01
      inHull1 = geometry::inhulln(datC1, as.matrix(toTest))
      toUse = profileLogLik[allPoints1,c('lognugget','profile')][inHull1,]
      
      interp1 = mgcv::gam(profile ~ s(lognugget, k=nrow(toUse), fx=TRUE), data=toUse)
      prof1 = data.frame(lognugget=seq(min(toUse[,1])-0.1, max(toUse[,1])+0.1, len=501))
      prof1$z = predict(interp1, prof1)
      
      plot(profileLogLik$lognugget, profileLogLik$profile, cex=.2, xlab="log(nugget)", ylab="profileLogL")
      points(toTest, col='red', cex=0.6)
      points(toUse, col='blue', cex=0.6, pch=3)
      lines(prof1$lognugget, prof1$z, col='green')
      abline(h =0, lty = 2, col='red')

   }
   
   
   
   
   if('anisoRatio' %in% paramToEstimate){
      result = as.data.table(cbind(LogLik, params[,"anisoRatio"]))
      colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 3) ,sep = ''), "anisoRatio")
      profileLogLik <- result[, .(profile=max(.SD)), by=anisoRatio]
      profileLogLik[,'profile'] <- profileLogLik[,'profile'] - breaks
      
      datC1= geometry::convhulln(profileLogLik)
      allPoints1 = unique(as.vector(datC1))
      toTest = profileLogLik[allPoints1,]
      toTest[,'profile'] = toTest[,'profile'] - 0.01
      inHull1 = geometry::inhulln(datC1, as.matrix(toTest))
      toUse = profileLogLik[allPoints1,][inHull1,]
      
      interp1 = mgcv::gam(profile ~ s(anisoRatio, k=nrow(toUse), fx=TRUE), data=toUse)
      prof1 = data.frame(anisoRatio=seq(min(toUse[,1])-0.1, max(toUse[,1])+0.1, len=501))
      prof1$z = predict(interp1, prof1)
      
      plot(profileLogLik$anisoRatio, profileLogLik$profile, cex=.2, xlab="anisoRatio", ylab="profileLogL")
      points(toTest, col='red', cex=0.6)
      points(toUse, col='blue', cex=0.6, pch=3)
      lines(prof1$anisoRatio, prof1$z, col='green')
      prof1 <- prof1[prof1$anisoRatio>0,]
      f1 <- approxfun(prof1$anisoRatio, prof1$z)
      #f1 <- approxfun(toUse$anisoRatio, toUse$profile)
      #curve(f1(x), add = TRUE, col = 'green', n = 1001)   #the number of x values at which to evaluate
      abline(h =0, lty = 2, col='red')
      lower = min(profileLogLik$anisoRatio)
      upper = max(profileLogLik$anisoRatio)
      MLE <- optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum
      
      ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
      abline(v =ci[1], lty = 2, col='red')
      abline(v =ci[2], lty = 2, col='red')
      if(length(ci)==1){
         if(ci > MLE){
            ci <- c(lower, ci)
            message("did not find lower ci for anisoRatio")
         }else{
            ci <- c(ci, upper)
            message("did not find upper ci for anisoRatio")}
      }
      
      if(length(ci)==0 | length(ci)>2){
         warning("error in params")
         ci <- c(NA, NA)
      }
      Table["anisoRatio",] <- c(MLE, ci)
      ####################### log plot ##############################
      # profileLogLik$loganisoRatio <- log(profileLogLik$anisoRatio)
      # datC1= geometry::convhulln(profileLogLik[,c('loganisoRatio','profile')])
      # allPoints1 = unique(as.vector(datC1))
      # toTest = profileLogLik[allPoints1,c('loganisoRatio','profile')]
      # toTest[,'profile'] = toTest[,'profile'] - 0.01
      # inHull1 = geometry::inhulln(datC1, as.matrix(toTest))
      # toUse = profileLogLik[allPoints1,c('loganisoRatio','profile')][inHull1,]
      # 
      # interp1 = mgcv::gam(profile ~ s(loganisoRatio, k=nrow(toUse), fx=TRUE), data=toUse)
      # prof1 = data.frame(loganisoRatio=seq(min(toUse[,1])-0.1, max(toUse[,1])+0.1, len=501))
      # prof1$z = predict(interp1, prof1)
      # 
      # plot(profileLogLik$loganisoRatio, profileLogLik$profile, cex=.2, xlab="log(anisoRatio)", ylab="profileLogL")
      # points(toTest, col='red', cex=0.6)
      # points(toUse, col='blue', cex=0.6, pch=3)
      # lines(prof1$loganisoRatio, prof1$z, col='green')
      # abline(h =0, lty = 2, col='red')
      
   }
   
   
   
   if('anisoAngleRadians' %in% paramToEstimate){
      result = as.data.table(cbind(LogLik, params[,"anisoAngleRadians"]))
      colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 3) ,sep = ''), "anisoAngleRadians")
      profileLogLik <- result[, .(profile=max(.SD)), by=anisoAngleRadians]
      profileLogLik[,'profile'] <- profileLogLik[,'profile'] - breaks
      
      datC1= geometry::convhulln(profileLogLik)
      allPoints1 = unique(as.vector(datC1))
      toTest = profileLogLik[allPoints1,]
      toTest[,'profile'] = toTest[,'profile'] - 18
      inHull1 = geometry::inhulln(datC1, as.matrix(toTest))
      toUse = profileLogLik[allPoints1,][inHull1,]
      toTest[,'profile'] = toTest[,'profile'] + 18
      
      interp1 = mgcv::gam(profile ~ s(anisoAngleRadians, k=nrow(toUse), fx=TRUE), data=toUse)
      prof1 = data.frame(anisoAngleRadians=seq(min(toUse[,1])-0.1, max(toUse[,1])+0.1, len=501))
      prof1$z = predict(interp1, prof1)
      
      plot(profileLogLik$anisoAngleRadians, profileLogLik$profile, cex=.2, xlab="anisoAngleRadians", ylab="profileLogL")
      points(toTest, col='red', cex=0.6)
      points(toUse, col='blue', cex=0.6, pch=3)
      lines(prof1$anisoAngleRadians, prof1$z, col='green')
      f1 <- approxfun(prof1$anisoAngleRadians, prof1$z)
      #f1 <- approxfun(toUse$anisoAngleRadians, toUse$profile)
      #curve(f1(x), add = TRUE, col = 'green', n = 1001)   #the number of x values at which to evaluate
      abline(h =0, lty = 2, col='red')
      lower = min(profileLogLik$anisoAngleRadians)
      upper = max(profileLogLik$anisoAngleRadians)
      MLE <- optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum
      
      ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
      abline(v =ci[1], lty = 2, col='red')
      abline(v =ci[2], lty = 2, col='red')
      if(length(ci)==1){
         if(ci > MLE){
            ci <- c(lower, ci)
            message("did not find lower ci for anisoAngleRadians")
         }else{
            ci <- c(ci, upper)
            message("did not find upper ci for anisoAngleRadians")}
      }
      
      if(length(ci)==0 | length(ci)>2){
         warning("error in params")
         ci <- c(NA, NA)
      }
      Table["anisoAngleRadians",] <- c(MLE, ci)
      ####################### log plot ##############################
      # profileLogLik$loganisoAngleRadians <- log(profileLogLik$anisoAngleRadians)
      # datC1= geometry::convhulln(profileLogLik[,c('loganisoAngleRadians','profile')])
      # allPoints1 = unique(as.vector(datC1))
      # toTest = profileLogLik[allPoints1,c('loganisoAngleRadians','profile')]
      # toTest[,'profile'] = toTest[,'profile'] - 0.01
      # inHull1 = geometry::inhulln(datC1, as.matrix(toTest))
      # toUse = profileLogLik[allPoints1,c('loganisoAngleRadians','profile')][inHull1,]
      # 
      # interp1 = mgcv::gam(profile ~ s(loganisoAngleRadians, k=nrow(toUse), fx=TRUE), data=toUse)
      # prof1 = data.frame(loganisoAngleRadians=seq(min(toUse[,1])-0.1, max(toUse[,1])+0.1, len=501))
      # prof1$z = predict(interp1, prof1)
      # 
      # plot(profileLogLik$loganisoAngleRadians, profileLogLik$profile, cex=.2, xlab="log(anisoAngleRadians)", ylab="profileLogL")
      # points(toTest, col='red', cex=0.6)
      # points(toUse, col='blue', cex=0.6, pch=3)
      # lines(prof1$loganisoAngleRadians, prof1$z, col='green')
      # abline(h =0, lty = 2, col='red')

   }
   
   
   
   
   
   
   
   if('anisoAngleDegrees' %in% paramToEstimate){
      result = as.data.table(cbind(LogLik, params[,"anisoAngleDegrees"]))
      colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 3) ,sep = ''), "anisoAngleDegrees")
      profileLogLik <- result[, .(profile=max(.SD)), by=anisoAngleDegrees]
      profileLogLik[,'profile'] <- profileLogLik[,'profile'] - breaks
      
      datC1= geometry::convhulln(profileLogLik)
      allPoints1 = unique(as.vector(datC1))
      toTest = profileLogLik[allPoints1,]
      toTest[,'profile'] = toTest[,'profile'] - 0.1
      inHull1 = geometry::inhulln(datC1, as.matrix(toTest))
      toUse = profileLogLik[allPoints1,][inHull1,]
      toTest[,'profile'] = toTest[,'profile'] + 0.1
      
      interp1 = mgcv::gam(profile ~ s(anisoAngleDegrees, k=nrow(toUse), fx=TRUE), data=toUse)
      prof1 = data.frame(anisoAngleDegrees=seq(min(toUse[,1])-0.1, max(toUse[,1])+0.1, len=501))
      prof1$z = predict(interp1, prof1)
      
      plot(profileLogLik$anisoAngleDegrees, profileLogLik$profile, cex=.2, xlab="anisoAngleDegrees", ylab="profileLogL")
      points(toTest, col='red', cex=0.6)
      points(toUse, col='blue', cex=0.6, pch=3)
      lines(prof1$anisoAngleDegrees, prof1$z, col='green')
      f1 <- approxfun(prof1$anisoAngleDegrees, prof1$z)
      #f1 <- approxfun(toUse$anisoAngleDegrees, toUse$profile)
      #curve(f1(x), add = TRUE, col = 'green', n = 1001)   #the number of x values at which to evaluate
      abline(h =0, lty = 2, col='red')
      lower = min(profileLogLik$anisoAngleDegrees)
      upper = max(profileLogLik$anisoAngleDegrees)
      MLE <- optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum
      
      ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
      abline(v =ci[1], lty = 2, col='red')
      abline(v =ci[2], lty = 2, col='red')
      if(length(ci)==1){
         if(ci > MLE){
            ci <- c(lower, ci)
            message("did not find lower ci for anisoAngleDegrees")
         }else{
            ci <- c(ci, upper)
            message("did not find upper ci for anisoAngleDegrees")}
      }
      
      if(length(ci)==0 | length(ci)>2){
         warning("error in params")
         ci <- c(NA, NA)
      }
      Table["anisoAngleDegrees",] <- c(MLE, ci)
      ####################### log plot ##############################
      # profileLogLik$loganisoAngleDegrees <- log(profileLogLik$anisoAngleDegrees)
      # datC1= geometry::convhulln(profileLogLik[,c('loganisoAngleDegrees','profile')])
      # allPoints1 = unique(as.vector(datC1))
      # toTest = profileLogLik[allPoints1,c('loganisoAngleDegrees','profile')]
      # toTest[,'profile'] = toTest[,'profile'] - 0.01
      # inHull1 = geometry::inhulln(datC1, as.matrix(toTest))
      # toUse = profileLogLik[allPoints1,c('loganisoAngleDegrees','profile')][inHull1,]
      # 
      # interp1 = mgcv::gam(profile ~ s(loganisoAngleDegrees, k=nrow(toUse), fx=TRUE), data=toUse)
      # prof1 = data.frame(loganisoAngleDegrees=seq(min(toUse[,1])-0.1, max(toUse[,1])+0.1, len=501))
      # prof1$z = predict(interp1, prof1)
      # 
      # plot(profileLogLik$loganisoAngleDegrees, profileLogLik$profile, cex=.2, xlab="log(anisoAngleDegrees)", ylab="profileLogL")
      # points(toTest, col='red', cex=0.6)
      # points(toUse, col='blue', cex=0.6, pch=3)
      # lines(prof1$loganisoAngleDegrees, prof1$z, col='green')
      # abline(h =0, lty = 2, col='red')  
     
     
   }
   
   
   
   
   ###############lambda hat#####################
   if(('boxcox'%in% paramToEstimate)  & length(boxcox)>5 ){
     likForboxcox = cbind(boxcox, apply(LogLik, 2,  max) )
     f1 <- approxfun(likForboxcox[,1], likForboxcox[,2]-breaks)
     # f1 <- splinefun(likForboxcox[,1], likForboxcox[,2]-breaks, method = "monoH.FC")
     plot(likForboxcox[,1], likForboxcox[,2]-breaks, ylab= "proLogL", xlab='boxcox', cex=0.5)
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

 
 
 # 2d
 result = as.data.table(cbind(LogLik, params[,c("anisoRatio","anisoAngleRadians")]))
 colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 3) ,sep = ''), "anisoRatio","anisoAngleRadians")
 profileLogLik <- result[, .(profile=max(.SD)), by=.(anisoRatio,anisoAngleRadians)]
 profileLogLik[,'profile'] <- profileLogLik[,'profile'] - breaks
 profileLogLik <- as.data.frame(profileLogLik)
 
 datC2 = geometry::convhulln(profileLogLik)
 allPoints2 = unique(as.vector(datC2))
 toTest2 = profileLogLik[allPoints2,]
 toTest2[,'profile'] = toTest2[,'profile'] - 0.1
 inHull2 = geometry::inhulln(datC2, as.matrix(toTest2))
 toUse2 = profileLogLik[allPoints2,][inHull2,]
 toTest2[,'profile'] = toTest2[,'profile'] + 0.1
 
 profileLogLik <- profileLogLik[order(profileLogLik[,1],decreasing=FALSE),]
 colDat2 = mapmisc::colourScale(profileLogLik[,'profile'], style='equal', dec=1, breaks=11, col='Spectral', rev=TRUE)
 plot(profileLogLik[,c('anisoRatio','anisoAngleRadians')], col = colDat2$plot, pch=16, cex=0.5)
 mapmisc::legendBreaks('topright', colDat2)
 points(toTest2[,c('anisoRatio','anisoAngleRadians')], col='black', cex=0.8)
 points(toUse2[,c('anisoRatio','anisoAngleRadians')], col='blue', pch=3, cex=0.8)
 

 interp2 = mgcv::gam(profile ~ s(anisoRatio,anisoAngleRadians, k=nrow(toUse2), fx=TRUE), data=toUse2)
 prof2list = list(anisoRatio=seq(min(toUse2[,1])-0.1, max(toUse2[,1])+0.1, len=101),
                  anisoAngleRadians=seq(min(toUse2[,2])-0.1, max(toUse2[,2])+0.1, len=101))
 prof2 = do.call(expand.grid, prof2list)
 prof2$z = predict(interp2, prof2)
 
 col2 = mapmisc::colourScale(prof2[,'z'], breaks= 11, dec=1, col='Spectral', rev=TRUE, style='equal')
 colPoints = mapmisc::colourScale(toUse2[,'profile'], breaks=col2$breaks, col=col2$col, style='fixed')
 plot(toUse2[,c('anisoRatio','anisoAngleRadians')], col=colPoints$plot, pch=15, cex=1)
 
 
 plot(prof2[,c('anisoRatio','anisoAngleRadians')], col=col2$plot, cex=2, pch=15)
 points(toUse2[,c('anisoRatio','anisoAngleRadians')], col='black')
 
 prof2 <- as.data.table(prof2)
 profileLogLik <- prof2[, .(profile=max(z)), by=anisoRatio]
 plot(profileLogLik[,c('anisoRatio','profile')],  cex=1, pch=15)
 plot(profileLogLik$anisoRatio, profileLogLik$profile, col='green')
 
 
 
 
 
 
 
 
 
 
 
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
             
             
             

             
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 