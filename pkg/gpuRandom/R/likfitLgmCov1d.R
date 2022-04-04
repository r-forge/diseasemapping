#' @title Estimate profile Log-likelihood for covariance parameters and lambda
#' @import data.table
#' @importFrom rootSolve uniroot.all
#' @useDynLib gpuRandom



get1dCovexhull <- function(profileLogLik,     # a data frame or data.table # 2 column names must be x1 and profile
                                a=0.1,    # minus a little thing
                                b=0,
                                m=1,
                                seqvalue){
  
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
  prof1 = data.frame(x1=seq(seqvalue[1], seqvalue[2], len=201))
  prof1$z = predict(interp1, prof1)
  
  output <- list(toUse = toUse, toTest = toTest, prof=prof1)
  
  output
}








#' @export
# betahat and sigmahat 
# profile log-likelihood for each covariance parameters + lambda
likfitLgmCov1d <- function(data,
                           formula, 
                           coordinates,
                           params, # CPU matrix for now, users need to provide proper parameters given their specific need
                           paramToEstimate = c('range','nugget'),
                           boxcox,  # boxcox is always estimated
                           cilevel,  # decimal
                           # seqNewRange,
                           # seqRange,
                           # seqShapelog,
                           # seqNugget,
                           # seqsdNugget,
                           # seqGamma3,
                           # seqGamma4,
                           # seqRatio,
                           # seqRadians,
                           type = c("float", "double")[1+gpuInfo()$double_support],
                           reml=FALSE, 
                           NparamPerIter,
                           Nglobal,
                           Nlocal,
                           NlocalCache,
                           verbose=FALSE){
  

  
  if(1 %in% boxcox){
    HasOne=TRUE
  }else{
    HasOne=FALSE 
  }
  
  
  if(0 %in% boxcox){
    boxcox = c(1, 0, setdiff(boxcox, c(1,0)))
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
    #XVYXVX3 <- as.matrix(XVYXVX)
    #which(is.na(XVYXVX2),arr.ind = TRUE)[,1]
    #which(is.na(XVYXVX3),arr.ind = TRUE)[,1]
    #tempp <- unique(which(is.na(XVYXVX2),arr.ind = TRUE)[,1])
    #tempp[-seq(1,length(tempp),2)]/2
    #unique(which(is.na(XVYXVX2),arr.ind = TRUE)[,1])
    
    a <- 0   
    for (j in 1:length(selected_rows)){
      a<-c(a, c(((selected_rows[j]-1)*Ncov+1): (selected_rows[j]*Ncov)))
    }
    a <- a[-1]
    XVYXVX2 <- XVYXVX2[-a,   ]
  }
  
  if(HasOne==FALSE){
    LogLikcpu <- as.matrix(LogLikcpu[,-1])
    ssqY2 <- ssqY2[,-1]
    ssqBetahat2 <- ssqBetahat2[,-1]
    ssqResidual2 <- ssqResidual2[,-1]
    Ndata <- length(boxcox) - 1
    XVYXVX2 <- XVYXVX2[,-1]
    jacobian <- as.vector(jacobian)[-1]
    boxcox <- boxcox[-1]
  }
  
  ############## output matrix ####################
  Table <- matrix(NA, nrow=length(paramToEstimate) + Ncov + 1, ncol=3)
  rownames(Table) <-  c(colnames(covariates), "sdSpatial", paramToEstimate)
  colnames(Table) <-  c("estimate", "lci", "uci")
  

  index <- which(LogLikcpu == max(LogLikcpu, na.rm = TRUE), arr.ind = TRUE)
  #################sigma hat#########################
  if(reml==FALSE)  {
    Table["sdSpatial",1] <- sqrt(ssqResidual2[index[1],index[2]]/Nobs)
  }else{         
    Table["sdSpatial",1] <- sqrt(ssqResidual2[index[1],index[2]]/(Nobs - Ncov))
  }
  
  
  maximum <- max(LogLikcpu)
  breaks = maximum - qchisq(cilevel,  df = 1)/2
  breaks2d = maximum - qchisq(cilevel,  df = 2)/2
  par(mfrow = c(3, 2))
  
  
  ############### profile for covariance parameters #####################
  if(params[1,'anisoRatio'] <= 1){
    gamma3 <-  unname(sqrt(1/paramsRenew[,'anisoRatio']-1) * cos(2*(paramsRenew[,'anisoAngleRadians'] + pi/2)))
    gamma4 <-  unname(sqrt(1/paramsRenew[,'anisoRatio']-1) * sin(2*(paramsRenew[,'anisoAngleRadians'] + pi/2)))
  }else{
    gamma3 <-  unname(sqrt(paramsRenew[,'anisoRatio']-1) * cos(2*(paramsRenew[,'anisoAngleRadians'])))
    gamma4 <-  unname(sqrt(paramsRenew[,'anisoRatio']-1) * sin(2*(paramsRenew[,'anisoAngleRadians'])))
  }
  aniso <- cbind(gamma3, gamma4)
  Newrange <- log(paramsRenew[,'range']^2/paramsRenew[,'anisoRatio'])
  paramsRenew <- cbind(paramsRenew, Newrange, aniso)
  #names(Newrange) <- 'Newrange'

  
  # Nconfig = c(NA,NA,NA,608, 740)[length(paramToEstimate)]
  # alphas = c(0.001,0.01, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 0.99, 0.999)
  # Salphas = c(NA, rep(alphas, each=Nconfig))
  # length(Salphas) == nrow(LogLikcpu)
  # colAlpha = mapmisc::colourScale(Salphas, style='unique', breaks = length(alphas), col='Spectral', opacity = 0.7)
  # colAlpha$plot[is.na(colAlpha$plot)] = '#000000FF'
  #xx = tapply(toPredictNatural$z, toPredictNatural$range, max)
  #par(mar = c(3,3,0.1, 0.1))
  ######################range ########
  if('Newrange' %in% paramToEstimate){
    result = as.data.table(cbind(Newrange, LogLikcpu))
    #head(result)
    colnames(result) <- c("x1", paste(c('boxcox'), round(boxcox, digits = 3) ,sep = ''))
    profileLogLik <- result[, .(profile=max(.SD)), by=.(x1)]
    profileLogLik[,'profile'] <- profileLogLik[,'profile'] - breaks
    #profileLogLik$col = colAlpha$plot

    profileLogLik <- profileLogLik[profile > maximum- breaks -10]  #
    profileLogLik <- as.data.frame(profileLogLik)

    #plot(exp(0.5*profileLogLik$x1), profileLogLik$profile, cex=.4, xlab="Newrange",pch=16, ylab="profileLogL", log='x') #,col=profileLogLik$col)
    plot(profileLogLik$x1, profileLogLik$profile, cex=.4, xlab="Newrange",pch=16, ylab="profileLogL") #,col=profileLogLik$col)
    #mapmisc::legendBreaks('right', colAlpha, bty='n')
    datC2 = geometry::convhulln(profileLogLik)
    allPoints = unique(as.vector(datC2))
    toTest = profileLogLik[allPoints,]
    toTest[,'profile'] = toTest[,'profile'] + 0.1
    inHull = geometry::inhulln(datC2, as.matrix(toTest))
    toUse = profileLogLik[allPoints,][!inHull,]
    toTest = profileLogLik[allPoints,]

    interp1 = mgcv::gam(profile ~ s(x1, k=nrow(toUse),  m=1, fx=TRUE), data=toUse)
    profNewrange = data.frame(x1=seq(min(toUse$x1), max(toUse$x1), len=501))
    profNewrange$z = predict(interp1, profNewrange)

    # points(exp(0.5*toTest[,1]), toTest[,2], col='red', cex=0.6)
    # points(exp(0.5*toUse[,1]), toUse[,2], col='blue', cex=0.6, pch=3)
    # lines(exp(0.5*profNewrange$x1), profNewrange$z, col = 'green')
    points(toTest[,1], toTest[,2], col='red', cex=0.6)
    points(toUse[,1], toUse[,2], col='blue', cex=0.6, pch=3)
    lines(profNewrange$x1, profNewrange$z, col = 'green')
    abline(h =0, lty = 2, col='red')
    lower = min(profileLogLik$x1)
    upper = max(profileLogLik$x1)
    f1 <- approxfun(profNewrange$x1, profNewrange$z)
    MLE <- optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum
    ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
    #abline(v =exp(0.5*c(MLE,ci)), lty = 2, col='red')
    abline(v =c(MLE,ci), lty = 2, col='red')
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
    Table["Newrange",] <- c(MLE,ci)

    ####################### log plot ##############################
    # profileLogLik$logrange <- log(profileLogLik$x1)
    # newdata <- profileLogLik[,c('logrange','profile')]
    # colnames(newdata)[1]<-"x1"
    # interlog <- get1dCovexhull(newdata)
    
    # plot(profileLogLik$logrange, profileLogLik$profile, cex=.2, xlab="log(range)", ylab="profileLogL")
    # points(interlog$toTest, col='red', cex=0.6)
    # points(interlog$toUse, col='blue', cex=0.6, pch=3)
    # lines(interlog$prof$x1, interlog$prof$z, col='green')
    # abline(h =0, lty = 2, col='red')
    
  }
  
  
  if('range' %in% paramToEstimate){
    result = data.table::as.data.table(cbind(LogLikcpu, paramsRenew[,"range"]))
    colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 3) ,sep = ''), "x1")
    profileLogLik <- result[, .(profile=max(.SD)), by=x1]
    profileLogLik[,'profile'] <- profileLogLik[,'profile'] - breaks
    profileLogLik <- profileLogLik[profile > maximum- breaks -10]
    #
    #  inter <- get1dCovexhull(profileLogLik, seqvalue = seqRange)    # a data frame or data.table # 2 column names must be x1 and profile
    # # #
    #   toTest <- inter$toTest
    #   toUse <- inter$toUse
    #   profRange <- inter$prof
    
    datC2 = geometry::convhulln(profileLogLik)
    allPoints = unique(as.vector(datC2))
    toTest = profileLogLik[allPoints,]
    toTest[,'profile'] = toTest[,'profile'] + 0.1
    inHull = geometry::inhulln(datC2, as.matrix(toTest))
    toUse = profileLogLik[allPoints,][!inHull,]
    toTest = profileLogLik[allPoints,]
    #
    interp1 = mgcv::gam(profile ~ s(x1, k=nrow(toUse),  m=1, fx=TRUE), data=toUse)
    profrange = data.frame(x1=seq(min(toUse$x1), max(toUse$x1), len=501))
    profrange$z = predict(interp1, profrange)
    
    
    plot(profileLogLik$x1, profileLogLik$profile, cex=.2, xlab="range", ylab="profileLogL")
    points(toTest, col='red', cex=0.6)
    points(toUse, col='blue', cex=0.6, pch=3)
    lines(profrange$x1, profrange$z, col = 'green')
    abline(h =0, lty = 2, col='red')
    lower = min(profileLogLik$x1)
    upper = max(profileLogLik$x1)
    f1 <- approxfun(profrange$x1, profrange$z)
    MLE <- optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum
    ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
    abline(v =c(MLE,ci), lty = 2, col='red')
    
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
    # profileLogLik$logrange <- log(profileLogLik$x1)
    # newdata <- profileLogLik[,c('logrange','profile')]
    # colnames(newdata)[1]<-"x1"
    # interlog <- get1dCovexhull(newdata)
    
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
    profileLogLik <- profileLogLik[profile > maximum- breaks -10]
    profileLogLik <- as.data.frame(profileLogLik)
    # inter <- get1dCovexhull(profileLogLik, seqvalue = c(0.1, 60))     # a data frame or data.table # 2 column names must be x1 and profile
    # 
    # toTest <- inter$toTest
    # toUse <- inter$toUse
    # prof1 <- inter$prof
    # 
    # plot(profileLogLik$x1, profileLogLik$profile, cex=.2, xlab="shape", ylab="profileLogL", xlim=c(0,30))
    #plot(log(profileLogLik$x1), profileLogLik$profile, cex=.2, xlab="logshape", ylab="profileLogL")
    # # points(toTest, col='red', cex=0.6)
    # # points(toUse, col='blue', cex=0.6, pch=3)
    # # f1 <- approxfun(inter$prof$x1, inter$prof$z)
    # lines(inter$prof$x1, inter$prof$z, col = 'green')
    # curve(f1(x), add = TRUE, col = 'green', n = 1001) 
    # abline(h =0, lty = 2, col='red')
    ####################### log plot ##############################
    profileLogLik$logshape <- log(profileLogLik$x1)
    newdata <- profileLogLik[,c('logshape','profile')]
    colnames(newdata)[1]<-"x1"
    #interlog <- get1dCovexhull(newdata,a=0.01,seqvalue = seqShapelog, m=1)
    
    datC2 = geometry::convhulln(newdata)
    allPoints = unique(as.vector(datC2))
    toTest = newdata[allPoints,]
    toTest[,'profile'] = toTest[,'profile'] + 0.1
    inHull = geometry::inhulln(datC2, as.matrix(toTest))
    toUse = newdata[allPoints,][!inHull,]
    toTest = newdata[allPoints,]
    toUse <- toUse[order(toUse$x1),]
    #toUse <- head(toUse, - 1)
    
    plot(profileLogLik$x1, profileLogLik$profile, cex=.2, xlab="shape", ylab="profileLogL",  log='x')
    points(exp(toTest[,1]),toTest[,2], col='red', cex=0.6)
    points(exp(toUse[,1]), toUse[,2], col='blue', cex=0.6, pch=3)
    
    interp1 = mgcv::gam(profile ~ s(x1, k=nrow(toUse),  m=1, fx=TRUE), data=toUse)
    profShapeLog = data.frame(x1=seq(min(toUse$x1), max(toUse$x1)-0.02, len=501))
    profShapeLog$z = predict(interp1, profShapeLog)
  
    #plot(newdata$x1, newdata$profile, cex=.2, xlab="log(shape)", ylab="profileLogL")

    lines(exp(profShapeLog$x1), profShapeLog$z, col = 'green')
    abline(h =0, lty = 2, col='red') 
    f1 <- approxfun(profShapeLog$x1, profShapeLog$z)
    #curve(f1(x), add = TRUE, col = 'green', n = 1001) 

    lower = min(toUse$x1)
    upper = max(toUse$x1)-0.05 
    MLE <- optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum
    ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
    #abline(v =c(MLE,ci), lty = 2, col='red')
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
    Table["shape",] <- c(exp(MLE), exp(ci))
    
    
  }       
  
  
  
  ################sd nugget ##############     
  
  
  if('sdNugget' %in% paramToEstimate){
    paramsRenew <- cbind(paramsRenew, sqrt(paramsRenew[,"nugget"]) * Table["sdSpatial",1])
    colnames(paramsRenew)[ncol(paramsRenew)] <- 'sdNugget'
    
    result = data.table::as.data.table(cbind(LogLikcpu, paramsRenew[,"sdNugget"]))
    colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 3) ,sep = ''), "x1")
    profileLogLik <- result[, .(profile=max(.SD)), by=x1]
    profileLogLik[,'profile'] <- profileLogLik[,'profile'] - breaks
    profileLogLik <- profileLogLik[profile > maximum- breaks -10]
    
    inter <- get1dCovexhull(profileLogLik, seqvalue = seqsdNugget)     # a data frame or data.table # 2 column names must be x1 and profile
    
    toTest <- inter$toTest
    toUse <- inter$toUse
    profsdNugget <- inter$prof
    
    # plot(profileLogLik$x1, profileLogLik$profile, cex=.2, xlab="sdNugget", ylab="profileLogL")
    # points(toTest, col='red', cex=0.6)
    # points(toUse, col='blue', cex=0.6, pch=3)
    # abline(h =0, lty = 2, col='red')
    profsdNugget <- profsdNugget[profsdNugget$x1>0,]
    f1 <- approxfun(profsdNugget$x1, profsdNugget$z)
    curve(f1(x), add = TRUE, col = 'green', n = 1001)
    lower = min(profsdNugget$x1)
    upper = max(profsdNugget$x1)
    MLE <- optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum
    ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
    #abline(v =ci[1], lty = 2, col='red')
    #abline(v =ci[2], lty = 2, col='red')
    
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
    result = data.table::as.data.table(cbind(LogLikcpu, paramsRenew[,"nugget"]))
    colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 3) ,sep = ''), "x1")
    profileLogLik <- result[, .(profile=max(.SD)), by=x1]
    profileLogLik[,'profile'] <- profileLogLik[,'profile'] - breaks
    profileLogLik <- profileLogLik[profile > maximum- breaks -10]
    
    plot(profileLogLik$x1, profileLogLik$profile, cex=.2, xlab="nugget", ylab="profileLogL")
    # inter <- get1dCovexhull(profileLogLik, seqvalue = seqNugget, m=1, a=0.00001)     # a data frame or data.table # 2 column names must be x1 and profile
    # toTest <- inter$toTest
    # toUse <- inter$toUse
    # profNugget <- inter$prof
    
    datC2 = geometry::convhulln(profileLogLik)
    allPoints = unique(as.vector(datC2))
    toTest = profileLogLik[allPoints,]
    toTest[,'profile'] = toTest[,'profile'] + 0.1
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
    if(nrow(toUse)>2){
    
    interp1 = mgcv::gam(profile ~ s(x1, k=nrow(toUse), m=1, fx=TRUE), data=toUse)
    profNugget = data.frame(x1=seq(min(toUse$x1), max(toUse$x1), len=501))
    profNugget$z = predict(interp1, profNugget)
    
    points(toTest, col='red', cex=0.6)
    points(toUse, col='blue', cex=0.6, pch=3)
    #profNugget <- profNugget[profNugget$x1>0,]
    f1 <- approxfun(profNugget$x1, profNugget$z)
    #toUse <- toUse[order(toUse$x1),]
    lines(profNugget$x1, profNugget$z, col = 'green')
    #lines(toUse$x1, toUse$profile, col = 'green')
    # curve(f1(x), add = TRUE, col = 'green', n = 1001)   #the number of x values at which to evaluate
    abline(h =0, lty = 2, col='red')
    # abline(h= maximum - qchisq(c(0.8, 0.95), df=1)/2 - breaks, col='red')
    lower = min(profileLogLik$x1)
    upper = max(profileLogLik$x1)
    MLE <- optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum
    ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
    abline(v =c(MLE,ci), lty = 2,  col='red')
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
    }else{
      interp1 = mgcv::gam(profile ~ s(x1, k=nrow(toUse), m=1, fx=TRUE), data=toUse)
      profNugget = data.frame(x1=seq(min(toUse$x1), max(toUse$x1), len=501))
      profNugget$z = predict(interp1, profNugget)
      
      #plot(profileLogLik$x1, profileLogLik$profile, cex=.2, xlab="nugget", ylab="profileLogL", log='x')
      points(toTest, col='red', cex=0.6)
      points(toUse, col='blue', cex=0.6, pch=3)
      #profNugget <- profNugget[profNugget$x1>0,]
      f1 <- approxfun(toUse$x1, toUse$profile)
      #toUse <- toUse[order(toUse$x1),]
      #lines(profNugget$x1, profNugget$z, col = 'green')
      lines(toUse$x1, toUse$profile, col = 'green')
      # curve(f1(x), add = TRUE, col = 'green', n = 1001)   #the number of x values at which to evaluate
      abline(h =0, lty = 2, col='red')
      # abline(h= maximum - qchisq(c(0.8, 0.95), df=1)/2 - breaks, col='red')
      lower = min(profileLogLik$x1)
      upper = max(profileLogLik$x1)
      MLE <- optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum
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
    }
    Table["nugget",] <- c(MLE, ci)
    ####################### log plot ##############################
    # profileLogLik$lognugget <- log(profileLogLik$x1)
    # newdata <- profileLogLik[,c('lognugget','profile')]
    # colnames(newdata)[1]<-"x1"
    # interlog <- get1dCovexhull(newdata,a=0.01, seqvalue=seqNuggetlog)
    
    # plot(profileLogLik$lognugget, profileLogLik$profile, cex=.2, xlab="log(nugget)", ylab="profileLogL")
    # points(interlog$toTest, col='red', cex=0.6)
    # points(interlog$toUse, col='blue', cex=0.6, pch=3)
    # lines(interlog$prof$x1, interlog$prof$z, col='green')
    # abline(h =0, lty = 2, col='red') 
    
  }
  
  if('anisoRatio' %in% paramToEstimate){
    result = as.data.table(cbind(LogLikcpu, paramsRenew[,"anisoRatio"]))
    colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 3) ,sep = ''), "x1")
    profileLogLik <- result[, .(profile=max(.SD)), by=x1]
    profileLogLik[,'profile'] <- profileLogLik[,'profile'] - breaks
    profileLogLik <- profileLogLik[profile > maximum- breaks -10]
    #
    datC2 = geometry::convhulln(profileLogLik)
    allPoints = unique(as.vector(datC2))
    toTest = profileLogLik[allPoints,]
    toTest[,'profile'] = toTest[,'profile'] + 0.1
    inHull = geometry::inhulln(datC2, as.matrix(toTest))
    toUse = profileLogLik[allPoints,][!inHull,]
    toTest = profileLogLik[allPoints,]

    interp1 = mgcv::gam(profile ~ s(x1, k=nrow(toUse), m=1, fx=TRUE), data=toUse)
    profRatio = data.frame(x1=seq(min(toUse$x1), max(toUse$x1), len=501))
    profRatio$z = predict(interp1, profRatio)
    #
    plot(profileLogLik$x1, profileLogLik$profile, cex=.2, xlab="anisoRatio", ylab="profileLogL")
    points(toTest, col='red', cex=0.6)
    points(toUse, col='blue', cex=0.6, pch=3)
    #lines(profRatio$x1, profRatio$z, col='green')
    f1 <- approxfun(profRatio$x1, profRatio$z)
    curve(f1(x), add = TRUE, col = 'green', n = 1001)
    abline(h =0, lty = 2, col='red')
    lower = min(profileLogLik$x1)
    upper = max(profileLogLik$x1)
    MLE <- optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum
    ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
    abline(v =c(MLE,ci), lty = 2, col='red')
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
    Table["anisoRatio",] <- c(MLE, ci)
  }
  
  
  if('anisoAngleRadians' %in% paramToEstimate){
    result = as.data.table(cbind(LogLikcpu, paramsRenew[,"anisoAngleRadians"]))
    colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 3) ,sep = ''), "x1")
    profileLogLik <- result[, .(profile=max(.SD)), by=x1]
    profileLogLik[,'profile'] <- profileLogLik[,'profile'] - breaks
    profileLogLik <- profileLogLik[profile >  maximum- breaks -10]
    #
    datC2 = geometry::convhulln(profileLogLik)
    allPoints = unique(as.vector(datC2))
    toTest = profileLogLik[allPoints,]
    toTest[,'profile'] = toTest[,'profile'] + 0.1
    inHull = geometry::inhulln(datC2, as.matrix(toTest))
    toUse = profileLogLik[allPoints,][!inHull,]
    toTest = profileLogLik[allPoints,]

    interp1 = mgcv::gam(profile ~ s(x1, k=nrow(toUse), m=1, fx=TRUE), data=toUse)
    profRadians = data.frame(x1=seq(min(toUse$x1), max(toUse$x1), len=501))
    profRadians$z = predict(interp1, profRadians)
    #
    plot(profileLogLik$x1, profileLogLik$profile, cex=.2, xlab="anisoAngleRadians", ylab="profileLogL")
    points(toTest, col='red', cex=0.6)
    points(toUse, col='blue', cex=0.6, pch=3)
    #toUse <- toUse[order(toUse[,1]),]
    #lines(toUse$x1, toUse$profile, col='green')
    lines(profRadians$x1, profRadians$z, col='green')
    f1 <- approxfun(profRadians$x1, profRadians$z)
    curve(f1(x), add = TRUE, col = 'green', n = 1001)
    abline(h =0, lty = 2, col='red')
    lower = min(profileLogLik$x1)
    upper = max(profileLogLik$x1)
    MLE <- optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum
    ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
    abline(v =c(MLE,ci), lty = 2, col='red')
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
    Table["anisoAngleRadians",] <- c(MLE, ci)
    ####################### log plot ##############################
    
  }
  
  
  if('gamma3' %in% paramToEstimate){
    #  result = as.data.table(cbind(LogLikcpu, aniso))
    #  colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 3) ,sep = ''), 'gamma3','gamma4')
    #  profileLogLik <- result[, .(profile=max(.SD)), by=.(gamma3, gamma4)]
    #  profileLogLik[,'profile'] <- profileLogLik[,'profile'] - breaks
    #  profileLogLik <- profileLogLik[profile > -10]
    #  profileLogLik <- as.data.frame(profileLogLik)
    # 
    #  datC2 = geometry::convhulln(profileLogLik)
    #  allPoints = unique(as.vector(datC2))
    #  toTest = profileLogLik[allPoints,]
    #  toTest[,'profile'] = toTest[,'profile'] + 0.1
    #  inHull = geometry::inhulln(datC2, as.matrix(toTest))
    #  toUse = profileLogLik[allPoints,][!inHull,]
    #  toTest = profileLogLik[allPoints,]
    #  
    #  
    #  
    #  fit = mgcv::gam(profile ~ s(gamma3, gamma4, k=nrow(toUse), m=1,fx=TRUE), data=toUse)
    #  toPredict = list(gamma3=seq(seqGamma3[1], seqGamma3[2], len=501),
    #                   gamma4=seq(seqGamma4[1], seqGamma4[2], len=501))
    #  toPredict = do.call(expand.grid, toPredict)
    #  toPredict$z = predict(fit, toPredict)
    # 
    #  yy = tapply(toPredict$z, toPredict$gamma3, max)
    #  plot(as.numeric(names(yy)), yy)
    # 
    #  qq = tapply(toPredict$z, toPredict$gamma4, max)
    #  plot(as.numeric(names(qq)), qq)
    #  abline(h =0, lty = 2, col='red')
    # 
    # col2 = mapmisc::colourScale(as.vector(toPredict[,'z']), breaks=SbreaksC, col=colDat2$col, style='fixed')
    # colPoints = mapmisc::colourScale(toUse2[,'profile'], breaks=SbreaksC, col=col2$col, style='fixed')
    # plot(toPredict[,c('gamma3','gamma4')], col=col2$plot, pch=15)
    # points(toUse2[,c('gamma3','gamma4')], col='black', cex=0.8)
    # mapmisc::legendBreaks('bottomright', breaks = Sprob, col=colDat2$col, bty='n')
    
    # prof2list = list(anisoRatio=seq(seqRatio[1], seqRatio[2], len=501),
    #                  anisoAngleRadians=seq(-0.5,0.5, len=501))
    # prof2natural = do.call(expand.grid, prof2list)
    # 
    # prof2naturalC = sqrt(prof2natural[,'anisoRatio']-1)* cos(2*(prof2natural[,'anisoAngleRadians'] + pi/2))+
    #   1i* sqrt(prof2natural[,'anisoRatio']-1)* sin(2*(prof2natural[,'anisoAngleRadians'] + pi/2))
    # 
    # 
    # prof2naturalAsGamma = data.frame(gamma3=Re(prof2naturalC), gamma4=Im(prof2naturalC))
    # prof2natural$z = predict(fit, prof2naturalAsGamma)
    
    # Sprob = c(0, 0.1, 0.2, 0.5, 0.8, 0.95, 0.99, 0.999, 1)
    # Sbreaks =qchisq(Sprob, df=2)/2
    # Sbreaks = pmin(Sbreaks, 1000)
    # Sbreaks[1] = -10
    # SbreaksC = rev(max(LogLikcpu) - breaks-Sbreaks)
    # colDat2 = mapmisc::colourScale(profileLogLik[,'profile'], style='fixed',
    #                                breaks=SbreaksC, 
    #                                col='Spectral', rev=TRUE)
    # colNatural = mapmisc::colourScale(prof2natural$z, col= colDat2$col, breaks=SbreaksC, style='fixed')
    # plot(prof2natural[,'anisoRatio'],prof2natural[,'anisoAngleRadians'], col=colNatural$plot, pch=16, ylab="anisoAngleRadians", xlab='anisoRatio')
    # mapmisc::legendBreaks('bottomright', breaks = rev(Sprob), col=colDat2$col, bty='n')
    # points(naturalspace, cex=0.8)
    
    # xx = tapply(prof2natural$z, prof2natural$anisoRatio, max)
    # plot(as.numeric(names(xx)), xx)

    # xx = tapply(prof2natural$z, prof2natural$anisoAngleRadians, max)
    # plot(as.numeric(names(xx)), xx) 
  
    result = as.data.table(cbind(LogLikcpu, aniso[,"gamma3"]))
    colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 3) ,sep = ''), "x1")
    profileLogLik <- result[, .(profile=max(.SD)), by=x1]
    profileLogLik[,'profile'] <- profileLogLik[,'profile'] - breaks
    profileLogLik <- profileLogLik[profile > maximum- breaks -10]
    
    # inter <- get1dCovexhull(profileLogLik, seqvalue = seqGamma3, m=1)     # a data frame or data.table # 2 column names must be x1 and profile
    # toTest <- inter$toTest
    # toUse <- inter$toUse
    # profGamma3 <- inter$prof
    
    datC2 = geometry::convhulln(profileLogLik)
    allPoints = unique(as.vector(datC2))
    toTest = profileLogLik[allPoints,]
    toTest[,'profile'] = toTest[,'profile'] + 0.1
    inHull = geometry::inhulln(datC2, as.matrix(toTest))
    toUse = profileLogLik[allPoints,][!inHull,]
    toTest = profileLogLik[allPoints,]
    
    interp1 = mgcv::gam(profile ~ s(x1, k=nrow(toUse), m=1, fx=TRUE), data=toUse)
    profGamma3 = data.frame(x1=seq(min(toUse$x1), max(toUse$x1), len=501))
    profGamma3$z = predict(interp1, profGamma3)
    
    plot(profileLogLik$x1, profileLogLik$profile, cex=.2, xlab="gamma3", ylab="profileLogL")
    points(toTest, col='red', cex=0.6)
    points(toUse, col='blue', cex=0.6, pch=3)
    
    #profGamma3 <- profGamma3[profGamma3$x1>0,]
    f1 <- approxfun(profGamma3$x1, profGamma3$z)
    lines(profGamma3$x1, profGamma3$z, col = 'green')
    #curve(f1(x), add = TRUE, col = 'green', n = 1001)   #the number of x values at which to evaluate
    abline(h =0, lty = 2, col='red')
    lower = min(profGamma3$x1)
    upper = max(profGamma3$x1)
    MLE <- optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum
    ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
    abline(v =c(MLE,ci), lty = 2, col='red')
    if(length(ci)==1){
      if(ci > MLE){
        ci <- c(lower, ci)
        message("did not find lower ci for gamma3")
      }else{
        ci <- c(ci, upper)
        message("did not find upper ci for gamma3")}
    }
    
    if(length(ci)==0 | length(ci)>2){
      warning("error in params")
      ci <- c(NA, NA)
    }
    Table["gamma3",] <- c(MLE, ci)
    ####################### log plot ##############################
    
  }
  
  
  
  if('gamma4' %in% paramToEstimate){
    result = as.data.table(cbind(LogLikcpu, aniso[,"gamma4"]))
    colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 3) ,sep = ''), "x1")
    profileLogLik <- result[, .(profile=max(.SD)), by=x1]
    profileLogLik[,'profile'] <- profileLogLik[,'profile'] - breaks
    profileLogLik <- profileLogLik[profile > maximum- breaks -10]
    
    plot(profileLogLik$x1, profileLogLik$profile, cex=.2, xlab="gamma4", ylab="profileLogL")
    
    datC2 = geometry::convhulln(profileLogLik)
    allPoints = unique(as.vector(datC2))
    toTest = profileLogLik[allPoints,]
    toTest[,'profile'] = toTest[,'profile'] + 0.1
    inHull = geometry::inhulln(datC2, as.matrix(toTest))
    toUse = profileLogLik[allPoints,][!inHull,]
    toTest = profileLogLik[allPoints,]
    
    interp1 = mgcv::gam(profile ~ s(x1, k=nrow(toUse), m=1, fx=TRUE), data=toUse)
    profGamma4 = data.frame(x1=seq(min(toUse$x1), max(toUse$x1), len=501))
    profGamma4$z = predict(interp1, profGamma4)
    
    # inter <- get1dCovexhull(profileLogLik, seqvalue = seqGamma4, a=0.00, m=1)     # a data frame or data.table # 2 column names must be x1 and profile
    # toTest <- inter$toTest
    # toUse <- inter$toUse
    # profGamma4 <- inter$prof
    
    points(toTest, col='red', cex=0.6)
    points(toUse, col='blue', cex=0.6, pch=3)
    
    #profGamma4 <- profGamma4[profGamma4$x1>0,]
    f1 <- approxfun(profGamma4$x1, profGamma4$z)
    lines(profGamma4$x1, profGamma4$z, col = 'green')
    #curve(f1(x), add = TRUE, col = 'green', n = 1001)   #the number of x values at which to evaluate
    abline(h =0, lty = 2, col='red')
    lower = min(profileLogLik$x1)
    upper = max(profileLogLik$x1)
    MLE <- optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum
    ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
    abline(v =c(MLE,ci), lty = 2, col='red')
    if(length(ci)==1){
      if(ci > MLE){
        ci <- c(lower, ci)
        message("did not find lower ci for gamma4")
      }else{
        ci <- c(ci, upper)
        message("did not find upper ci for gamma4")}
    }
    
    if(length(ci)==0 | length(ci)>2){
      warning("error in params")
      ci <- c(NA, NA)
    }
    Table["gamma4",] <- c(MLE, ci)
    ####################### log plot ##############################
    
  }
  
  # if('anisoAngleDegrees' %in% paramToEstimate){
  #   result = as.data.table(cbind(LogLikcpu, paramsRenew[,"anisoAngleDegrees"]))
  #   colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 3) ,sep = ''), "x1")
  #   profileLogLik <- result[, .(profile=max(.SD)), by=x1]
  #   profileLogLik[,'profile'] <- profileLogLik[,'profile'] - breaks
  #   inter <- get1dCovexhull(profileLogLik)     # a data frame or data.table # 2 column names must be x1 and profile
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
    plot(likForboxcox[,1], likForboxcox[,2]-breaks, ylab= "proLogL", xlab='boxcox', cex=0.5)
    curve(f1(x), add = TRUE, col = 2, n = 1001)   #the number of x values at which to evaluate
    abline(h =0, lty = 2)
    
    lower = min(boxcox)
    upper = max(boxcox)
    MLE <- optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum
    ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
    abline(v =c(MLE,ci), lty = 2)

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
  }
  
  NcolTotal = Ndata + Ncov
  
  ###############betahat#####################
  Betahat <- matrix(0, nrow=Ncov, ncol=Ndata)
  a<-c( ((index[1]-1)*Ncov+1) : (index[1]*Ncov) )
  mat <- XVYXVX2[a,((Ndata+1):NcolTotal)]
  mat[upper.tri(mat)] <- mat[lower.tri(mat)]
  Betahat <- solve(mat) %*% XVYXVX2[a,index[2]]
  Table[colnames(covariates), 1] <- Betahat
  
  # Betahat <- rep(0, Ncov)
  # for(i in 1:Nparam){
  #   a<-c( ((i-1)*Ncov+1) : (i*Ncov) )
  #   mat <- XVYXVX2[a,((Ndata+1):NcolTotal)]
  #   mat[upper.tri(mat)] <- mat[lower.tri(mat)]
  #   Betahat0 <- t(solve(mat) %*% XVYXVX2[a,index[2]])
  #   Betahat <- rbind(Betahat,Betahat0)
  # }
  # # 
  #  x<- cbind(Betahat[-1,], LogLikcpu)
  #  colnames(x)[1:Ncov] <- c('intercept',colnames(covariates)[-1])
  #  x <- as.data.frame(x)
   # xSub = x[x$boxcox1 > (maximum - 10), ]
   
   #####
   # selected_rows2 <- which(x$boxcox1 > (maximum - 10))
   # a <- 0   
   # for (j in 1:length(selected_rows2)){
   #   a<-c(a, c(((selected_rows2[j]-1)*Ncov+1): (selected_rows2[j]*Ncov)))
   # }
   # a <- a[-1]
   # XTVinvX3 <- XVYXVX2[a, (ncol(XVYXVX2)-Ncov+1):ncol(XVYXVX2)  ]
   # #####
   # head( XTVinvX3)
   # appendBetaMatrix <- matrix(0, nrow=2*nrow(xSub), ncol=4)
   # colnames(appendBetaMatrix) <- colnames(x)
   # 
   # 
   # intercept <- seq(0.5,3.5, len=100)
   # quadratic = 0.5*(- covarianceInv[1,1]  ) * (intercept-xSub[i,'intercept'])^2 + maximum
   # lines(intercept, quadratic, col="green")
   # points(xSub[i,'intercept'] + 2*interceptsd,  0.5*(- covarianceInv[1,1] ) * (2*interceptsd)^2 + maximum, col='red')
   # 
   # 
   # for(i in 1:nrow(xSub)){
   #   interval <- c(((i-1)*Ncov+1) : (i*Ncov))
   #   covarianceInv <- XTVinvX[interval, ]
   #   covariance <- solve(covarianceInv)
   #   interceptsd <- sqrt(covariance[1,1])
   #   cov1sd <- sqrt(covariance[2,2])
   #   #onesd <- xSub[i, 'intercept'] + interceptsd
   #   
   #   LogLikI <- 0.5*(-covarianceInv[1,1]) * (2*interceptsd)^2 + maximum
   #   appendBetaMatrix[2*i-1, 1:2] <- c( unlist(xSub[i,'intercept']) + 2*interceptsd, LogLikI)
   #   appendBetaMatrix[2*i, 1:2] <- c( unlist(xSub[i,'intercept']) - 2*interceptsd, LogLikI)
   #   
   #   LogLikC <- 0.5*(-covarianceInv[2,2]) * (2*cov1sd)^2 + maximum
   #   appendBetaMatrix[2*i-1, 3:4] <- c( unlist(xSub[i,'cov1']) + 2*cov1sd, LogLikC)
   #   appendBetaMatrix[2*i, 3:4] <- c( unlist(xSub[i,'cov1']) - 2*cov1sd, LogLikC)
   # 
   #   #qnorm(c(0.09, 0.999), mean = x[i, 'intercept'], sd = interceptsd)
   #   #LogLik_i <- 0.5* 2*c(interceptsd,cov1sd) %*% (-covarianceInv) %*% (2*c(interceptsd,cov1sd)) + maximum
   #   #appendBetaMatrix[2*i-1,] <- c( unlist(xSub[i, c('intercept', 'cov1')]) + 2*c(interceptsd,  cov1sd), LogLik_i)
   #   #appendBetaMatrix[2*i,] <- c( unlist(xSub[i, c('intercept', 'cov1')]) - 2*c(interceptsd,  cov1sd), LogLik_i)
   # }
   # head(appendBetaMatrix)
   # colnames(appendBetaMatrix) <- c('intercept','boxcox1' ,'cov1', 'boxcox1')
   # head(xSub)
   # 
   # 
   # xfinalI <- rbind(xSub[,c(1,3)], appendBetaMatrix[,1:2])
   # plot(xfinal$cov1, xfinal$boxcox1, cex=0.2)
   # plot(xfinalI$intercept, xfinalI$boxcox1, cex=0.2)
   
  #  plot(xSub$cov1, xSub$boxcox1, cex=0.2)
  #  plot(xSub$intercept, xSub$boxcox1 - breaks, cex=0.2)
  #  abline(h = 0, lty = 2, col='red')
  #  abline(v = c(0.4660031, 3.196358), lty = 2, col='red')
  #  abline(v = c(0.8043255,  2.849402), col='purple')
  # # plot(xSub$cov2, xSub$boxcox1, cex=0.2)
  #  plot(x$`(Intercept)`, x$boxcox1, cex=0.2) 
  

  
  if(all(c('sdNugget','shape', 'gamma3') %in% paramToEstimate)){
  Output <- list(LogLik=LogLikcpu,
                 breaks = breaks,
                 breaks2d = breaks2d,
                 mleIndex = index,
                 summary = Table,
                 #BetahatTable = x,
                 profNewrange = profNewrange,
                 profShapeLog = profShapeLog,
                 profNugget = profNugget,
                 profsdNugget = profsdNugget,
                 profGamma3 = profGamma3,
                 profGamma4 = profGamma4,
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
  }
  
  if(!('shape' %in% paramToEstimate)){
    Output <- list(LogLik=LogLikcpu,
                   breaks = breaks,
                   breaks2d = breaks2d,
                   mleIndex = index,
                   summary = Table,
                   #BetahatTable = x,
                   profNewrange = profNewrange,
                   profNugget = profNugget,
                   profGamma3 = profGamma3,
                   profGamma4 = profGamma4,
                   params = paramsRenew,
                   Infindex = selected_rows,
                   boxcox = boxcox,
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
  
    }
  
  
 if('anisoRatio' %in% paramToEstimate & ('anisoAngleRadians' %in% paramToEstimate)){
      Output <- list(LogLik=LogLikcpu,
                   breaks = breaks,
                   breaks2d = breaks2d,
                   mleIndex = index,
                   summary = Table,
                   #BetahatTable = x,
                   profrange = profrange,
                   profShapeLog = profShapeLog,
                   profNugget = profNugget,
                   profRatio = profRatio,
                   profRadians = profRadians,
                   params = paramsRenew,
                   Infindex = selected_rows,
                   boxcox = boxcox,
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
  }
  
  if('gamma3' %in% paramToEstimate & ('gamma4' %in% paramToEstimate)){
    Output <- list(LogLik=LogLikcpu,
                   breaks = breaks,
                   breaks2d = breaks2d,
                   mleIndex = index,
                   summary = Table,
                   #BetahatTable = x,
                   profNewrange = profNewrange,
                   profShapeLog = profShapeLog,
                   profNugget = profNugget,
                   profGamma3 = profGamma3,
                   profGamma4 = profGamma4,
                   params = paramsRenew,
                   Infindex = selected_rows,
                   boxcox = boxcox,
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
  }  
  Output
  
  
  
}


