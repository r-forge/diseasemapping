#' @title Estimate profile Log-likelihood for covariance parameters and lambda
#' @import data.table
#' @importFrom rootSolve uniroot.all
#' @useDynLib gpuRandom



get2dCovexhullinter <- function(profileLogLik,     # a data frame or data.table # 2 column names must be x1 and profile
                                a=0.1,    # minus a little thing
                                seqvalue,
                                m=1){
  
  datC2 = geometry::convhulln(profileLogLik)
  allPoints = unique(as.vector(datC2))
  toTest = profileLogLik[allPoints,]
  toTest[,'profile'] = toTest[,'profile'] + a
  inHull = geometry::inhulln(datC2, as.matrix(toTest))
  toUse = profileLogLik[allPoints,][!inHull,]
  toTest = profileLogLik[allPoints,]
  
  # datC1= geometry::convhulln(profileLogLik)
  # allPoints1 = unique(as.vector(datC1))
  # toTest = profileLogLik[allPoints1,]
  # toTest[,'profile'] = toTest[,'profile'] - a
  # toTest[,'x1'] = toTest[,'x1'] + b
  # inHull1 = geometry::inhulln(datC1, as.matrix(toTest))
  # toUse = profileLogLik[allPoints1,][inHull1,]
  # toTest[,'profile'] = toTest[,'profile'] + a
  # toTest[,'x1'] = toTest[,'x1'] - b
  
  interp1 = mgcv::gam(profile ~ s(x1, k=nrow(toUse), m=m, fx=TRUE), data=toUse)
  prof1 = data.frame(x1=seq(seqvalue[1], seqvalue[2], len=101))
  prof1$z = predict(interp1, prof1)
  
  output <- list(toUse = toUse, toTest = toTest, prof=prof1)
  
  output
}


#' @export
# betahat and sigmahat 
# profile log-likelihood for each covariance parameters + lambda
 likfitLgmCov2d <- function(data,
                         formula, 
                         coordinates,
                         params, # CPU matrix for now, users need to provide proper parameters given their specific need
                         paramToEstimate = c('range','nugget'),
                         boxcox,  # boxcox is always estimated
                         cilevel,  # decimal
                         df = 2,
                         seqNewRange,
                         seqShapelog,
                         seqNugget,
                         seqRatio,
                         seqRadians,
                         type = c("float", "double")[1+gpuInfo()$double_support],
                         reml=FALSE, 
                         NparamPerIter,
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

  # params0[which(is.na(as.vector(detVar))),]

  
  if(reml== FALSE){ 
    # if(fixVariance == FALSE){ # ml
    gpuRandom:::matrix_vector_sumBackend(Nobs*log(ssqResidual/Nobs),
                                         detVar,
                                         jacobian,  
                                         Nobs + Nobs*log(2*pi),
                                         minusTwoLogLik,
                                         Nglobal)

  }else if(reml==TRUE){
    # if(fixVariance == FALSE){# remlpro
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
  
  selected_rows <- which(is.na(as.vector(detVar)))
  if(length(selected_rows)==0){
    paramsRenew <- params
    detVar2 <- as.vector(detVar)
    detReml2 <- as.vector(detReml)
    ssqY2 <- as.matrix(ssqY)
    ssqBetahat2 = as.matrix(ssqBetahat)
    ssqResidual2 = as.matrix(ssqResidual)
    XVYXVX2 <- as.matrix(XVYXVX)
  }else{
    Nparam = Nparam - length(selected_rows)
    paramsRenew <- params[-selected_rows,]
    LogLikcpu <- as.matrix(LogLikcpu[-selected_rows,]) 
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
  
  ############## output matrix ####################
  Table <- matrix(NA, nrow=length(paramToEstimate) + Ncov + 1, ncol=1)
  rownames(Table) <-  c(colnames(covariates), "sdSpatial", paramToEstimate)
  colnames(Table) <-  c("estimate")
  
  
  index <- which(LogLikcpu == max(LogLikcpu, na.rm = TRUE), arr.ind = TRUE)
  #################sigma hat#########################
  if(reml==FALSE)  {
    Table["sdSpatial",1] <- sqrt(ssqResidual[index[1],index[2]]/Nobs)
  }else{         
    Table["sdSpatial",1] <- sqrt(ssqResidual[index[1],index[2]]/(Nobs - Ncov))
  }
  
  
  maximum <- max(LogLikcpu)
  breaks = maximum - qchisq(cilevel,  df = 1)/2
  breaksoutput = maximum - qchisq(cilevel,  df = df)/2
  par(mfrow = c(2, 2))
  
  
  ############### profile for covariance parameters #####################
  if(params[1,'anisoRatio'] <= 1){
    gamma3 <-  unname(sqrt(1/paramsRenew['anisoRatio']-1) * cos(2*(paramsRenew['anisoAngleRadians'] + pi/2)))
    gamma4 <-  unname(sqrt(1/paramsRenew['anisoRatio']-1) * sin(2*(paramsRenew['anisoAngleRadians'] + pi/2)))
  }else{
    gamma3 <-  unname(sqrt(paramsRenew[,'anisoRatio']-1) * cos(2*(paramsRenew[,'anisoAngleRadians'])))
    gamma4 <-  unname(sqrt(paramsRenew[,'anisoRatio']-1) * sin(2*(paramsRenew[,'anisoAngleRadians'])))
  }
  aniso <- cbind(gamma3, gamma4)
  
  Newrange <- log(paramsRenew[,'range']^2/paramsRenew[,'anisoRatio'])
  paramsRenew <- cbind(paramsRenew, Newrange, aniso)
  ############### profile for covariance parameters #####################
  ######################range ########
  if('Newrange' %in% paramToEstimate){
    result = as.data.table(cbind(Newrange, LogLikcpu))
    #head(result)
    colnames(result) <- c("x1", paste(c('boxcox'), round(boxcox, digits = 3) ,sep = ''))
    profileLogLik <- result[, .(profile=max(.SD)), by=.(x1)]
    profileLogLik[,'profile'] <- profileLogLik[,'profile'] - breaks
    profileLogLik <- profileLogLik[profile > maximum- breaks -10]
    profileLogLik <- as.data.frame(profileLogLik)
    
    plot(exp(0.5*profileLogLik$x1), profileLogLik$profile, cex=.2, xlab="Newrange", ylab="profileLogL", log='x')
    
    datC2 = geometry::convhulln(profileLogLik)
    allPoints = unique(as.vector(datC2))
    toTest = profileLogLik[allPoints,]
    toTest[,'profile'] = toTest[,'profile'] + 0.1
    inHull = geometry::inhulln(datC2, as.matrix(toTest))
    toUse = profileLogLik[allPoints,][!inHull,]
    toTest = profileLogLik[allPoints,]
    
    interp1 = mgcv::gam(profile ~ s(x1, k=nrow(toUse),  m=1, fx=TRUE), data=toUse)
    profNewrange = data.frame(x1=seq(seqNewRange[1], seqNewRange[2], len=101))
    profNewrange$z = predict(interp1, profNewrange)
    
    points(exp(0.5*toTest[,1]), toTest[,2], col='red', cex=0.6)
    points(exp(0.5*toUse[,1]), toUse[,2], col='blue', cex=0.6, pch=3)
    lines(exp(0.5*profNewrange$x1), profNewrange$z, col = 'green')
    abline(h =0, lty = 2, col='red')
    lower = min(profileLogLik$x1)
    upper = max(profileLogLik$x1)
    f1 <- approxfun(profNewrange$x1, profNewrange$z)
    MLE <- optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum
    ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
    abline(v =exp(0.5*c(MLE,ci)), lty = 2, col='red') 
    if(length(ci)==1){
      if(ci > MLE){
        ci <- c(lower, ci)
        message("did not find lower ci for Newrange")
      }else{
        ci <- c(ci, upper)
        message("did not find upper ci for Newrange")}
    }
    
    if(length(ci)==0 | length(ci)>2){
      warning("error in params")
      ci <- c(NA, NA)
    }
    Table["Newrange",] <- c(MLE)
    # result = data.table::as.data.table(cbind(LogLikcpu, paramsRenew[,"range"]))
    # colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 3) ,sep = ''), "x1")
    # profileLogLik <- result[, .(profile=max(.SD)), by=x1]
    # profileLogLik[,'profile'] <- profileLogLik[,'profile'] - breaks
    # profileLogLik <- profileLogLik[profile > -20]
    # 
    # inter <- get2dCovexhullinter(profileLogLik)    # a data frame or data.table # 2 column names must be x1 and profile
    # 
    # toTest <- inter$toTest
    # toUse <- inter$toUse
    # prof1 <- inter$prof
    # 
    # # plot(profileLogLik$x1, profileLogLik$profile, cex=.2, xlab="range", ylab="profileLogL")
    # # points(toTest, col='red', cex=0.6)
    # # points(toUse, col='blue', cex=0.6, pch=3)
    # f1 <- approxfun(prof1$x1, prof1$z)
    # # curve(f1(x), add = TRUE, col = 'green', n = 1001)
    # # abline(h =0, lty = 2, col='red')
    # lower = min(profileLogLik$x1)
    # upper = max(profileLogLik$x1)
    # MLE <- optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum
    # Table["range",] <- MLE
    ####################### log plot ##############################
    # profileLogLik$logrange <- log(profileLogLik$x1)
    # newdata <- profileLogLik[,c('logrange','profile')]
    # colnames(newdata)[1]<-"x1"
    # interlog <- get2dCovexhullinter(newdata)
    
    # plot(profileLogLik$logrange, profileLogLik$profile, cex=.2, xlab="log(range)", ylab="profileLogL")
    # points(interlog$toTest, col='red', cex=0.6)
    # points(interlog$toUse, col='blue', cex=0.6, pch=3)
    # lines(interlog$prof$x1, interlog$prof$z, col='green')
    # abline(h =0, lty = 2, col='red')

  }
  
  
  
  
  ################shape ##############   
  if('shape' %in% paramToEstimate){
    result = data.table::as.data.table(cbind(LogLikcpu, paramsRenew[,"shape"]))
    colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 3) ,sep = ''), "x1")
    profileLogLik <- result[, .(profile=max(.SD)), by=x1]
    profileLogLik[,'profile'] <- profileLogLik[,'profile'] - breaks
    profileLogLik <- profileLogLik[profile > maximum - breaks -10]
    profileLogLik <- as.data.frame(profileLogLik)
    plot(profileLogLik$x1, profileLogLik$profile, cex=.2, xlab="shape", ylab="profileLogL",  log='x')
    profileLogLik$logshape <- log(profileLogLik$x1)
    newdata <- profileLogLik[,c('logshape','profile')]
    colnames(newdata)[1]<-"x1"
    interlog <- get2dCovexhullinter(newdata,a=0.1,seqvalue = seqShapelog, m=1)
    profShapeLog <- interlog$prof
    points(exp(interlog$toTest[,1]),interlog$toTest[,2], col='red', cex=0.6)
    points(exp(interlog$toUse[,1]), interlog$toUse[,2], col='blue', cex=0.6, pch=3)
    lines(exp(profShapeLog$x1), profShapeLog$z, col = 'green')
    f1 <- approxfun(profShapeLog$x1, profShapeLog$z)
    #curve(f1(x), add = TRUE, col = 'green', n = 1001) 
    abline(h =0, lty = 2, col='red') 
    lower = log(min(profileLogLik$x1))
    upper =8
    MLE <- optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum
    ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
    abline(v =exp(c(MLE,ci)), lty = 2, col='red')
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
    Table["shape",] <- c(exp(MLE))
    # result = data.table::as.data.table(cbind(LogLikcpu, paramsRenew[,"shape"]))
    # colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 3) ,sep = ''), "x1")
    # profileLogLik <- result[, .(profile=max(.SD)), by=x1]
    # profileLogLik[,'profile'] <- profileLogLik[,'profile'] - breaks
    # profileLogLik <- profileLogLik[profile > -20]
    # 
    # # inter <- get2dCovexhullinter(profileLogLik)     # a data frame or data.table # 2 column names must be x1 and profile
    # # 
    # # toTest <- inter$toTest
    # # toUse <- inter$toUse
    # # prof1 <- inter$prof
    # 
    # # plot(profileLogLik$x1, profileLogLik$profile, cex=.2, xlab="shape", ylab="profileLogL", xlim=c(0,100))
    # # points(toTest, col='red', cex=0.6)
    # # points(toUse, col='blue', cex=0.6, pch=3)
    # # f1 <- approxfun(inter$prof$x1, inter$prof$z)
    # # curve(f1(x), add = TRUE, col = 'green', n = 1001) 
    # # abline(h =0, lty = 2, col='red')
    # ####################### log plot ##############################
    # profileLogLik$logshape <- log(profileLogLik$x1)
    # newdata <- profileLogLik[,c('logshape','profile')]
    # colnames(newdata)[1]<-"x1"
    # interlog <- get2dCovexhullinter(newdata,a=0.2)
    # 
    # # plot(newdata$x1, newdata$profile, cex=.2, xlab="log(shape)", ylab="profileLogL")
    # # points(interlog$toTest, col='red', cex=0.6)
    # # points(interlog$toUse, col='blue', cex=0.6, pch=3)
    # f1 <- approxfun(interlog$prof$x1, interlog$prof$z)
    # # curve(f1(x), add = TRUE, col = 'green', n = 1001) 
    # # abline(h =0, lty = 2, col='red') 
    # lower = log(min(profileLogLik$x1))
    # upper = log(50)
    # MLE <- optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum
    # Table["shape",] <- c(exp(MLE))
    
    
  }       
  
  
  
  ################sd nugget ##############     
 if('sdNugget' %in% paramToEstimate){
    paramsRenew <- cbind(paramsRenew, sqrt(paramsRenew[,"nugget"]) * Table["sdSpatial",1])
    colnames(paramsRenew)[ncol(paramsRenew)] <- 'sdNugget'
    
    result = data.table::as.data.table(cbind(LogLikcpu, paramsRenew[,"sdNugget"]))
    colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 3) ,sep = ''), "x1")
    profileLogLik <- result[, .(profile=max(.SD)), by=x1]
    profileLogLik[,'profile'] <- profileLogLik[,'profile'] - breaks
    profileLogLik <- profileLogLik[profile > -20]
    
    inter <- get2dCovexhullinter(profileLogLik)     # a data frame or data.table # 2 column names must be x1 and profile
    
    toTest <- inter$toTest
    toUse <- inter$toUse
    prof1 <- inter$prof
    
    # plot(profileLogLik$x1, profileLogLik$profile, cex=.2, xlab="sdNugget", ylab="profileLogL")
    # points(toTest, col='red', cex=0.6)
    # points(toUse, col='blue', cex=0.6, pch=3)
    # abline(h =0, lty = 2, col='red')
    prof1 <- prof1[prof1$x1>0,]
    f1 <- approxfun(prof1$x1, prof1$z)
    # curve(f1(x), add = TRUE, col = 'green', n = 1001)
    lower = min(prof1$x1)
    upper = max(prof1$x1)
    MLE <- optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum
    

    Table["sdNugget",] <- c(MLE)
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
    result = data.table::as.data.table(cbind(LogLikcpu, paramsRenew[,"nugget"]))
    colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 3) ,sep = ''), "x1")
    profileLogLik <- result[, .(profile=max(.SD)), by=x1]
    profileLogLik[,'profile'] <- profileLogLik[,'profile'] - breaks
    profileLogLik <- profileLogLik[profile > maximum- breaks-10]
    
    datC2 = geometry::convhulln(profileLogLik)
    allPoints = unique(as.vector(datC2))
    toTest = profileLogLik[allPoints,]
    toTest[,'profile'] = toTest[,'profile'] + 0.1
    inHull = geometry::inhulln(datC2, as.matrix(toTest))
    toUse = profileLogLik[allPoints,][!inHull,]
    toTest = profileLogLik[allPoints,]
    
    interp1 = mgcv::gam(profile ~ s(x1, k=nrow(toUse), m=1, fx=TRUE), data=toUse)
    profNugget = data.frame(x1=seq(seqNugget[1], seqNugget[2], len=101))
    profNugget$z = predict(interp1, profNugget)
 
    plot(profileLogLik$x1, profileLogLik$profile, cex=.2, xlab="nugget", ylab="profileLogL")
    points(toTest, col='red', cex=0.6)
    points(toUse, col='blue', cex=0.6, pch=3)
    f1 <- approxfun(profNugget$x1, profNugget$z)
    lines(profNugget$x1, profNugget$z, col = 'green')
    abline(h =0, lty = 2, col='red')
    lower = min(profileLogLik$x1)
    upper = max(profileLogLik$x1)
    MLE <- optimize(f1, c(lower, 3), maximum = TRUE, tol = 0.0001)$maximum
    ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
    abline(v =c(0,ci), lty = 2,  col='red')
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
    Table["nugget",] <- c(MLE)
    # result = data.table::as.data.table(cbind(LogLikcpu, paramsRenew[,"nugget"]))
    # colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 3) ,sep = ''), "x1")
    # profileLogLik <- result[, .(profile=max(.SD)), by=x1]
    # profileLogLik[,'profile'] <- profileLogLik[,'profile'] - breaks
    # profileLogLik <- profileLogLik[profile > -20]
    # 
    # inter <- get2dCovexhullinter(profileLogLik)     # a data frame or data.table # 2 column names must be x1 and profile
    # toTest <- inter$toTest
    # toUse <- inter$toUse
    # prof1 <- inter$prof
    # 
    # # plot(profileLogLik$x1, profileLogLik$profile, cex=.2, xlab="nugget", ylab="profileLogL")
    # # points(toTest, col='red', cex=0.6)
    # # points(toUse, col='blue', cex=0.6, pch=3)
    # prof1 <- prof1[prof1$x1>0,]
    # f1 <- approxfun(prof1$x1, prof1$z)
    # # curve(f1(x), add = TRUE, col = 'green', n = 1001)   #the number of x values at which to evaluate
    # # abline(h =0, lty = 2, col='red')
    # lower = min(prof1$x1)
    # upper = max(prof1$x1)
    # MLE <- optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum
    # # ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
    # # abline(v =ci[1], lty = 2, col='red')
    # # abline(v =ci[2], lty = 2, col='red')
    # # if(length(ci)==1){
    # #   if(ci > MLE){
    # #     ci <- c(lower, ci)
    # #     message("did not find lower ci for nugget")
    # #   }else{
    # #     ci <- c(ci, upper)
    # #     message("did not find upper ci for nugget")}
    # # }
    # # 
    # # if(length(ci)==0 | length(ci)>2){
    # #   warning("error in params")
    # #   ci <- c(NA, NA)
    # # }
    # Table["nugget",] <- c(MLE)
    ####################### log plot ##############################
    # profileLogLik$lognugget <- log(profileLogLik$x1)
    # newdata <- profileLogLik[,c('lognugget','profile')]
    # colnames(newdata)[1]<-"x1"
    # interlog <- get2dCovexhullinter(newdata,a=0.01)
    
    # plot(profileLogLik$lognugget, profileLogLik$profile, cex=.2, xlab="log(nugget)", ylab="profileLogL")
    # points(interlog$toTest, col='red', cex=0.6)
    # points(interlog$toUse, col='blue', cex=0.6, pch=3)
    # lines(interlog$prof$x1, interlog$prof$z, col='green')
    # abline(h =0, lty = 2, col='red') 
    
  }
  
  
  
  
  if(('gamma3' %in% paramToEstimate) | ('anisoRatio' %in% paramToEstimate)){
    result = as.data.table(cbind(LogLikcpu, paramsRenew[,"anisoRatio"]))
    colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 3) ,sep = ''), "x1")
    profileLogLik <- result[, .(profile=max(.SD)), by=x1]
    profileLogLik[,'profile'] <- profileLogLik[,'profile'] - breaks
    profileLogLik <- profileLogLik[profile > maximum- breaks-10]
    
    inter <- get2dCovexhullinter(profileLogLik, seqvalue = seqRatio)     # a data frame or data.table # 2 column names must be x1 and profile
    toTest <- inter$toTest
    toUse <- inter$toUse
    prof1 <- inter$prof
    
    plot(profileLogLik$x1, profileLogLik$profile, cex=.2, xlab="anisoRatio", ylab="profileLogL")
    points(toTest, col='red', cex=0.6)
    points(toUse, col='blue', cex=0.6, pch=3)

    f1 <- approxfun(prof1$x1, prof1$z)
     curve(f1(x), add = TRUE, col = 'green', n = 1001)   #the number of x values at which to evaluate
     abline(h =0, lty = 2, col='red')
    lower = min(prof1$x1)
    upper = max(prof1$x1)
    MLE <- optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum
    Table["anisoRatio",] <- c(MLE)
    ####################### log plot ##############################
    
  }
  
  
  
  if('anisoAngleRadians' %in% paramToEstimate){
    result = as.data.table(cbind(LogLikcpu, paramsRenew[,"anisoAngleRadians"]))
    colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 3) ,sep = ''), "x1")
    profileLogLik <- result[, .(profile=max(.SD)), by=x1]
    profileLogLik[,'profile'] <- profileLogLik[,'profile'] - breaks
    profileLogLik <- profileLogLik[profile > maximum- breaks-5]
    
    inter <- get2dCovexhullinter(profileLogLik, seqvalue = seqRadians)     # a data frame or data.table # 2 column names must be x1 and profile
    toTest <- inter$toTest
    toUse <- inter$toUse
    prof1 <- inter$prof
    
    plot(profileLogLik$x1, profileLogLik$profile, cex=.2, xlab="anisoAngleRadians", ylab="profileLogL")
    points(toTest, col='red', cex=0.6)
    points(toUse, col='blue', cex=0.6, pch=3)
    lines(prof1$x1, prof1$z, col='green')
    f1 <- approxfun(prof1$x1, prof1$z)
    abline(h =0, lty = 2, col='red')
    lower = min(profileLogLik$x1)
    upper = max(profileLogLik$x1)
    MLE <- optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum
    Table["anisoAngleRadians",] <- c(MLE)
    ####################### log plot ##############################
    
  }
  
  
  
  
  
  
  
  # if('anisoAngleDegrees' %in% paramToEstimate){
  #   result = as.data.table(cbind(LogLikcpu, paramsRenew[,"anisoAngleDegrees"]))
  #   colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 3) ,sep = ''), "x1")
  #   profileLogLik <- result[, .(profile=max(.SD)), by=x1]
  #   profileLogLik[,'profile'] <- profileLogLik[,'profile'] - breaks
  #   profileLogLik <- profileLogLik[profile > -20]
  #   
  #   inter <- get2dCovexhullinter(profileLogLik)     # a data frame or data.table # 2 column names must be x1 and profile
  #   toTest <- inter$toTest
  #   toUse <- inter$toUse
  #   prof1 <- inter$prof
  #   
  #   # plot(profileLogLik$x1, profileLogLik$profile, cex=.2, xlab="anisoAngleDegrees", ylab="profileLogL")
  #   # points(toTest, col='red', cex=0.6)
  #   # points(toUse, col='blue', cex=0.6, pch=3)
  #   # lines(prof1$x1, prof1$z, col='green')
  #   f1 <- approxfun(prof1$x1, prof1$z)
  #   # abline(h =0, lty = 2, col='red')
  #   lower = min(profileLogLik$x1)
  #   upper = max(profileLogLik$x1)
  #   MLE <- optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum
  #   Table["anisoAngleDegrees",] <- c(MLE)
  # 
  #   
  #   
  # }
  
  
  
  
  ###############lambda hat#####################
  if(('boxcox'%in% paramToEstimate)  & length(boxcox)>5 ){
    likForboxcox = cbind(boxcox, apply(LogLikcpu, 2,  max) )
    f1 <- approxfun(likForboxcox[,1], likForboxcox[,2]-breaks)
    # plot(likForboxcox[,1], likForboxcox[,2]-breaks, ylab= "proLogL", xlab='boxcox', cex=0.5)
    # curve(f1(x), add = TRUE, col = 2, n = 1001)   #the number of x values at which to evaluate
    # abline(h =0, lty = 2)

    lower = min(boxcox)
    upper = max(boxcox)
    MLE <- optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum

    Table["boxcox",] <- c(MLE)
  }else if(is.element('boxcox',paramToEstimate)  & length(boxcox) <= 5){
    message("boxcox: not enough values for interpolation!")
  }else if(is.element('boxcox',paramToEstimate)){
    Table["boxcox",1] <- boxcox[index[2]]
  }
  
  NcolTotal = Ndata + Ncov
  
  ###############betahat#####################
  Betahat <- matrix(0, nrow=Ncov, ncol=Ndata)
  a<-c( ((index[1]-1)*Ncov+1) : (index[1]*Ncov) )
  mat <- XVYXVX[a,((Ndata+1):NcolTotal)]
  mat[upper.tri(mat)] <- mat[lower.tri(mat)]
  Betahat <- solve(mat) %*% XVYXVX[a,index[2]]
  
  Table[colnames(covariates), 1] <- Betahat
  

  Output <- list(LogLik=LogLikcpu,
                 breaks = breaksoutput,
                 summary = Table,
                 params = paramsRenew,
                 Infindex = selected_rows,
                 Nobs = Nobs,
                 Ncov = Ncov,
                 Ndata = Ndata,
                 Nparam = Nparam,
                 ssqY = as.matrix(ssqY2),     
                 ssqBetahat = ssqBetahat2,
                 ssqResidual = ssqResidual2,
                 detVar = as.vector(detVar2),   
                 detReml = as.vector(detReml2),   
                 jacobian = as.vector(jacobian),
                 XVYXVX = as.matrix(XVYXVX2))
  
  
  
  Output
  
  
  
}













