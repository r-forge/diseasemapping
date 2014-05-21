objfunM = function(oparam,distVec,sqrtVar,
                   range=NULL,shape=NULL){
  
  if(!is.null(range))
    oparam['range'] = as.vector(range)
  if(!is.null(shape))
    oparam['shape'] = as.vector(shape)
  if(!any(names(oparam)=='variance')){	
    oparam['variance']=1
  }
  if(!any(names(oparam)=='nugget')){	
    oparam['nugget']=1
  }
  
  theM =
    geostatsp::matern(
      x=distVec, 
      param=oparam[c('range','shape','variance')]
    )
  thezeros = which(distVec<=0)
  theM[thezeros] = theM[thezeros] + oparam['nugget']
  
  sum( (sqrt(theM) - sqrtVar)^2)
}

loglikGmrfGivenQ = function(
  propNugget, 
  Yp, YL, 
  Xp, XL,
  Qchol, detQ, 
  QLL=NULL,cholQLL=NULL,
  XLp = NULL, YLp = NULL,
  empirical=NULL) {
  
  
  # propNugget =    tau^2/xi^2
  
  N= length(Yp) + c(ml=0,reml=-ncol(Xp))
  
  if(propNugget>0){
    
    
    # Q = L L', QLL = L'L
    # cholIcQ = chol of L'L + I xisq/tausq
    if(is.null(QLL)) {
      QLL = crossprod(expand(Qchol)$L)
      
      cholIcQ = Cholesky(QLL, LDL=FALSE,
                         Imult=1/propNugget)
      YLp = solve(cholIcQ, YL,
                  system='P')	
      XLp = solve(cholIcQ, XL,
                  system='P')	
    } else {
      cholIcQ = update(cholQLL, QLL,mult=1/propNugget)
    }
    
    
    Ybreve = solve(cholIcQ, 
                   YLp,	system='L')
    Xbreve = solve(cholIcQ, 
                   XLp, system='L')	
    
    XprecXinv = solve(Matrix::crossprod(Xbreve,Xbreve))
    
    betaHat = as.vector(
      XprecXinv %*% 
        Matrix::crossprod(Xbreve,Ybreve) 
    )
    names(betaHat) = colnames(XL)
    
    residsL =  as.vector(Ybreve -  Xbreve %*% betaHat )
    rpr = as.numeric(crossprod(residsL,residsL))
    
    
    
    logDetVar = as.numeric(
      2*determinant(cholIcQ,logarithm=TRUE)$modulus -
        detQ 
    )
    
    
    
    tausq = constHat= N/rpr
    xisq = tausq/propNugget
    
    m2logL = logDetVar + N*log(rpr) 
    
    
    
  } else { # no nugget 
    
    XprecX = Matrix::crossprod(XL,XL)
    XprecXinv = solve(XprecX)
    
    betaHat = as.numeric(XprecXinv %*% 
                           Matrix::crossprod(XL , YL))
    names(betaHat) = colnames(XL)
    
    residsL  = as.vector(YL - XL %*% betaHat)
    rpr = as.numeric(crossprod(residsL,residsL))
    
    logDetVar=-detQ
    tausq=0*N
    xisq = constHat = rpr/N
    m2logL = logDetVar + N*log(rpr)
    # matrix operations
    # Qchol, cholCovMat = Matrix::chol(covMat)
    # XL, cholCovInvX = Matrix::solve(cholCovMat, covariates)
    # XprecX, cholCovInvXcross = Matrix::crossprod(cholCovInvX)
    # XprecXinv, cholCovInvXcrossInv = Matrix::solve(cholCovInvXcross)
  }
  
  variances = c(tausq = tausq, xisq = xisq)
  
  
  m2logL['reml'] =	m2logL['reml'] +  determinant(
    XprecXinv,logarithm=TRUE)$modulus
  
  names(m2logL) = paste("m2logL.",names(m2logL),sep='')
  
  logL = -m2logL/2
  names(logL) = gsub("^m2","",names(logL))
  
  sebeta = diag(XprecXinv)
  sebeta = rep(sebeta,2)*
    rep(constHat,
        rep(length(sebeta),2))
  
  names(sebeta)= paste('se',rep(names(betaHat),2), 
                       names(sebeta),sep='.')
  
  varbetahat = abind::abind(
    as.matrix(XprecXinv)/constHat[1],
    as.matrix(XprecXinv)/constHat[2],
    along=3)
  
  dimnames(varbetahat) = 
    list(names(betaHat),names(betaHat),
         names(constHat))
  
  if(FALSE){#}!is.null(empirical)){
    
    empirical$yMult = empirical$yMult * variances['sigmasq.ml']	
    empirical[empirical$x<=0,'yMult'] = 
      empirical[empirical$x<=0,'yMult'] + 
      variances['tausq.ml']
    startparam = c(
      shape=2,
      range=2,
      variance=as.vector(variances['sigmasq.ml'])
    )
    if(variances['tausq.ml']>0)
      startparam = c(startparam,
                     nugget=as.vector(variances['tausq.ml'])				
      )
    
    postfit=optim(fn=objfunM,par=startparam,
                  lower=c(shape=0.1,range=0.1,
                          variance=0,nugget=0)[names(startparam)],
                  upper=c(shape=3,range=12,
                          variance=Inf,nugget=Inf)[names(startparam)],
                  distVec=empirical$x, 
                  sqrtVar = sqrt(empirical$yMult))$par
    
    
    
    names(postfit) = paste(names(postfit), ".postfit",sep='')
    
  } else {
    postfit=NULL
  }		
  
  result =c(m2logL, logL,
            variances,
            beta=betaHat, 
            sebeta,
            logDetVar=as.numeric(logDetVar),
            rpr=rpr,postfit)
  
  
  attributes(result)$varbetahat = varbetahat
  
  result
  
}

loglikGmrfOneRange = function(
  oneminusar,
  Yvec, Xmat, NN, propNugget=0,
  shape=1,  
  adjustEdges=FALSE,adjustParam=FALSE,
  adjustShape=FALSE,adjustMarginalVariance=FALSE) {
    
    NN =  maternGmrfPrec(NN,
                         param=c(shape=as.vector(shape),
                                 oneminusar=as.vector(oneminusar[1]),
                                 conditionalVariance=1),
                         adjustEdges=adjustEdges,adjustParam=adjustParam,
                         adjustShape=adjustShape,
                         adjustMarginalVariance=adjustMarginalVariance)
    
	
  
  cholPrec = Cholesky(NN,LDL=FALSE)
  LofQ = expand(cholPrec)$L
  
  theDet = 2*as.numeric(determinant(cholPrec,
                                    logarithm=TRUE)$modulus)
  
  Xp = solve(cholPrec,Xmat,system="P")
  Yp = solve(cholPrec,Yvec,system="P")
  
  
  
  XL = Matrix::crossprod(LofQ,Xp)
  YL = as.vector(
    Matrix::crossprod(LofQ,Yp)
  )
  
  argList = list(Yp=Yp,YL=YL,
                 Xp=Xp,XL=XL,
                 Qchol=cholPrec,
                 detQ=theDet,
                 empirical=attributes(NN)$param$empirical)
  
  if(!all(propNugget==0)) {
    # calculate t(LofQ) LofQ
    argList$QLL = crossprod(expand(argList$Qchol)$L)
    
    argList$cholQLL = Cholesky(argList$QLL, LDL=FALSE)
    argList$YLp =   solve(argList$cholQLL , YL,
                          system='P')	
    argList$XLp = solve(argList$cholQLL , XL,
                        system='P')	
    
  }
  
  if(FALSE) {
    # variance of Q inverse
    theraster = attributes(NN)$raster
    Nx=ncol(theraster);Ny=nrow(theraster)
    
    midcellCoord = c(round(Nx*.5),round(Ny*0.5)) # the middle cell
    midcell = c(Nx*(Ny-midcellCoord[2]) + midcellCoord[1])
    midVec = sparseMatrix(midcell,1,x=1,
                          dims=c(ncol(NN),1))
    
    Qcentre = solve(cholPrec,
                    midVec)[midcell]
    
  } else {
    Qcentre = NULL
  }
  
  if(length(propNugget)==1) {
    
    argList$propNugget=propNugget
    
    res = do.call(loglikGmrfGivenQ,argList)
    
    res = c(res,
            oneminusar=as.numeric(oneminusar),
            propNugget=as.numeric(propNugget),
            shape=as.numeric(shape),
            Qcentre=as.numeric(Qcentre)
    )
    
  } else {
    
    res = mapply(loglikGmrfGivenQ,
                 propNugget=propNugget,			
                 MoreArgs=argList,SIMPLIFY=TRUE
    )
    colnames(res) = paste("propNugget=", propNugget,sep="")
    
    res = rbind(res,
                oneminusar=as.numeric(oneminusar),
                propNugget=as.numeric(propNugget),
                shape=as.numeric(shape),
                Qcentre=as.numeric(Qcentre)
    )
    
  }
  
  attributes(res)$Qinfo = attributes(NN)$param
  
  res
  
}



loglikGmrf = function(
  oneminusar,
  Yvec, Xmat, NN, propNugget=0,
  shape=1, 
  adjustEdges=FALSE,adjustParam=FALSE,
  adjustShape=FALSE, 
  adjustMarginalVariance=FALSE,
  mc.cores=1) {
  
  
  argList = list(Yvec=Yvec, 
                 Xmat=Xmat, NN=NN,
                 propNugget=propNugget,
                 shape=shape,
                 adjustEdges=adjustEdges,
                 adjustParam=adjustParam,
                 adjustShape=adjustShape,
                 adjustMarginalVariance=adjustMarginalVariance
  )
  
  if(mc.cores>1) {
    myapply= function(...){
      parallel::mcmapply(...,
                         mc.cores=mc.cores)
    }
  } else {
    myapply= mapply
  }
  
    res=myapply(loglikGmrfOneRange,
                oneminusar=oneminusar,
                MoreArgs=argList,
                SIMPLIFY=FALSE
    )
    names(res) = 
      paste("oneminusar=",oneminusar,sep="")
  
  if(class(res[[1]])!= 'character'){
    Qinfo = attributes(res[[1]])$Qinfo
    
    res= simplify2array(res)
    
    attributes(res)$Qinfo = Qinfo
    
  }
  res
}




summaryGmrfFit= function(x, cellSize=1) {
  UseMethod("summaryGmrfFit")	
}

summaryGmrfFit.array = function(x, cellSize=1) {
  
  
  x2 = aperm(x,c(3,2,1))
  x2 = matrix(c(x2), ncol=dim(x2)[3])
  colnames(x2) = dimnames(x)[[1]]
  res = summaryGmrfFit.matrix(x2,npar=2)
  
  res$profL = list()
  res$profL$propNugget = cbind(propNugget=x['propNugget',,1],
                               logL.ml=apply(x['logL.ml',,],1,max),
                               logL.reml = apply(x['logL.reml',,],1,max)
  )
  whichL = apply(x['logL.ml',,],2,which.max)
  
  res$profL$oneminusar = cbind(
    oneminusar = x['oneminusar',1,],
    logL.ml=apply(x['logL.ml',,],2,max),
    logL.reml=apply(x['logL.reml',,],2,max)
  )
  
  if(all(
    c('range.postfit','shape.postfit') %in%
      dimnames(x)[[1]])) {
    res$profL$oneminusar  = cbind(
      res$profL$oneminusar,
      range.postit = diag(
        x['range.postfit',whichL
          ,]),
      shape.postit = diag(
        x['shape.postfit',whichL
          ,])
    )
  }
  
  
  #do something with the range and sigma squared
  # from oneminusar and xisq and shape
  rangge = cellSize*sqrt(2*x['shape',1,]*(1-x['oneminusar',1,])/x['oneminusar',1,])
    res$profL$range = list(
      ml = cbind(
        range =rangge, 
        logL=  apply(x['logL.ml',,],2,max)
        ),
      reml = cbind(
        range = rangge, 
      logL=  apply(x['logL.reml',,],2,max)
    )
      ) 
  
  
  res$profL$rangeInCells = list(
    ml = cbind(
      rangeInCells=	rangge/cellSize,
      logL=	apply(x['logL.ml',,],2,max)
    ),
    reml = cbind(
      rangeInCells=rangge/cellSize,
      logL=apply(x['logL.reml',,],2,max)
    )
  )
  

  
  res$profL$sigmasq = list(
    ml = cbind(
      sigmasq = x['xisq.ml',1,]/(4^(x['shape',1,]+1)*pi*x['shape',1,])
        *(1-x['oneminusar',1,])^x['shape',1,]/(x['oneminusar',1,])^x['shape',1,],
      logL=  apply(x['logL.ml',,],2,max)
      ),
    reml = cbind(
      sigmasq = x['xisq.reml',1,]/(4^(x['shape',1,]+1)*pi*x['shape',1,])
      *(1-x['oneminusar',1,])^x['shape',1,]/(x['oneminusar',1,])^x['shape',1,],
      logL=  apply(x['logL.reml',,],2,max)
    )
  )
  return(res)	
  
}

summaryGmrfFit.matrix = function(x,npar=1, cellSize=1) {
  
  if(any(rownames(x)=='logL.ml')){
    x = t(x)
  }
  result = list()
  someL = c('ml','reml')
  for(D1 in someL) {
    D=paste('logL.',D1,sep='')
    MLErow = x[which.max(x[ ,D]),]	
    thebetas = grep("^beta",names(MLErow), value=TRUE)
    betaMLE = MLErow[thebetas]	
    
    thebetas2 = gsub("^beta\\.","",thebetas)
    betaSE = MLErow[paste("se.",thebetas2,".ml",sep="")]
    
    betamat = cbind(mle=betaMLE,
                    se=betaSE,
                    q0.025=betaMLE - 1.96*betaSE,
                    q0.975=betaMLE + 1.96*betaSE,
                    pval = 2*pnorm(
                      abs(betaMLE)/
                        betaSE,
                      lower=FALSE)	
    )
    
    withinCI = which(-2*(x[,D] - 
                           MLErow[D]) < qchisq(0.95, 2))
    
    xsub = x[,grep("^se\\.",colnames(x),invert=TRUE),drop=FALSE]
    parCI= t(apply(
      xsub[withinCI, ,drop=FALSE],
      2,range))
    colnames(parCI) = c('q0.025', 'q0.975')
    
    parMat = cbind(mle=MLErow[colnames(xsub)],parCI,se=NA,pval=NA)
    parMat[rownames(betamat),colnames(betamat)] = betamat
    
    notD = someL[someL != D1]
    notD = grep(
      paste("\\.",notD, "$",sep=''),
      rownames(parMat), 
      invert=TRUE,value=TRUE)
    
    parMat = parMat[notD,]
    rownames(parMat) = gsub(paste("\\.", D1, "$",sep=""),
                            "",rownames(parMat))
    
    result[[D1]] = parMat
  }
  return(result)		
}