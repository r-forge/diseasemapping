

loglikGmrfGivenQ = function(
		propNugget,
  Ry, Y, 
  Rx, X,
  Q, Qchol=NULL, 
  detQ=NULL) {
  
  
  # propNugget =    tau^2/xi^2
  
  N= length(Y) + c(ml=0,reml=-ncol(X))

  if(is.null(Qchol)) {
	  Qchol = Cholesky(Q, LDL=FALSE)
  }
  if(is.null(detQ)) {
	  detQ= 2*determinant(Qchol,logarithm=TRUE)$modulus
  }
  
  if(propNugget>0){

	  xisqtausq = 1/propNugget
  
	  Vchol = update(Qchol, Q, mult=xisqtausq)
	  
	  Vy = solve(Vchol, Y)
	  Vx = solve(Vchol, X)
	  
    
    XprecXinv = solve(Matrix::crossprod(Rx,Vx))
    
    betaHat = as.vector(
      XprecXinv %*% 
        Matrix::crossprod(Rx,Vy) 
    )
    names(betaHat) = colnames(X)

	R = as.numeric(
			Matrix::crossprod(Ry,Vy) - 
					Matrix::crossprod(Ry, Vx) %*% betaHat
	)

 
	
	constHat=tausq = N/R
	xisq = tausq/propNugget
	
    logDetVar = as.numeric(
      2*determinant(Vchol,logarithm=TRUE)$modulus
	)
        
  } else { # no nugget 
    
	  
    XprecX = Matrix::crossprod(X,Rx)
    XprecXinv = solve(XprecX)
    
    betaHat = as.numeric(XprecXinv %*% 
                           Matrix::crossprod(X , Ry))
    names(betaHat) = colnames(X)
    
    resids  = as.vector(Y - X %*% betaHat)
    R = as.numeric(crossprod(resids,Q) %*% resids)
    
    logDetVar=0
    tausq=0
	
    constHat = xisq =R/N
    # matrix operations
    # Qchol, cholCovMat = Matrix::chol(covMat)
    # XL, cholCovInvX = Matrix::solve(cholCovMat, covariates)
    # XprecX, cholCovInvXcross = Matrix::crossprod(cholCovInvX)
    # XprecXinv, cholCovInvXcrossInv = Matrix::solve(cholCovInvXcross)
  }
  
  logDetVar = logDetVar - detQ + N - N*log(N)
  m2logL = logDetVar + N*log(R)
  
  
  variances = c(tausq = tausq, xisq = xisq,
		  sigmasq = xisq * attributes(Q)$param$theo['variance'])
  
  
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
  
   result =c(m2logL, logL,
            variances,
            beta=betaHat, 
            sebeta,
            logDetVar=as.numeric(logDetVar))
  
  
  attributes(result)$varbetahat = varbetahat
  
  result
  
}

loglikGmrfOneRange = function(
  oneminusar,
  Yvec, Xmat, NN, propNugget=0,
  shape=1,  
  adjustEdges=FALSE,adjustParam=FALSE) {
    
    Q =  maternGmrfPrec(NN,
                         param=c(shape=as.vector(shape),
                                 oneminusar=as.vector(oneminusar[1]),
                                 conditionalVariance=1),
		adjustEdges=adjustEdges,adjustParam=adjustParam)
    
	
  
	Qchol = Cholesky(Q,LDL=FALSE)
 
  
  detQ = 2*as.numeric(determinant(Qchol,
                                    logarithm=TRUE)$modulus)
  
  Rx = Q %*% Xmat
  Ry = Q %*% Yvec
  
  argList = list(
  	Ry=Ry, Y=Yvec, 
  	Rx=Rx, X=Xmat,
 	Q=Q, Qchol=Qchol, 
  detQ=detQ)  
  
thepar=attributes(Q)$paramInfo
   
  if(length(propNugget)==1) {
    
    argList$propNugget=propNugget
    
    res = do.call(loglikGmrfGivenQ,argList)
    
    res = c(res,
            oneminusar=as.numeric(oneminusar),
            propNugget=as.numeric(propNugget),
			shape = thepar$theo['shape'],
			shapeOpt = thepar$optimal['shape'],
			range = thepar$theo['range'],
			rangeOpt = thepar$optimal['range']
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
				shape = thepar$theo['shape'],
				shapeOpt = thepar$optimal['shape'],
				range = thepar$theo['range'],
				rangeOpt = thepar$optimal['range']
    )
    
  }
  
  attributes(res)$Qinfo = thepar
  
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
                 adjustParam=adjustParam
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
