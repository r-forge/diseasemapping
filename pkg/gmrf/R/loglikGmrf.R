bceps = 0.01

logLbc = function(bc, y, x, 
    logy=log(y), 
    xvx = crossprod(x),     
    xvxinv = solve(xvx),
    sumlogy = sum(logy)
){

  if(abs(bc)< bceps){
    y = logy
  } else if (abs(bc-1)>bceps){
     y = exp(bc*logy)-1
  }
  twologj = -2*(sumlogy)*(bc-1)

  
  betaHat = xvxinv %*% crossprod(x, y)
  ssq = y - x %*% betaHat
  sum(ssq^2) + twologj
  
}

loglikGmrfOneRange = function(
    oneminusar,
    Yvec, Xmat, NN, propNugget=0,
    shape=1,
    boxcox=1,
    fix.boxcox=TRUE,
    Nboxcox=5,
    reml=TRUE,
    sumLogY = NULL,
    adjustEdges=FALSE,
    optimizer=FALSE) {
  

  Q =  maternGmrfPrec(NN,
      param=c(shape=as.vector(shape),
          oneminusar=as.vector(oneminusar[1]),
          conditionalVariance=1),
      adjustEdges=adjustEdges)
  
  thepar=attributes(Q)$param
  
  propNugget = sort(unique(c(0, propNugget)))
  
  Ny = ncol(Yvec)
  if(is.null(Ny))
    Ny = 1
  
  if(!fix.boxcox){
    if(Ny != 1) warning('cant do box-cox with more than one dataset')
    if(any(Yvec)<0) warning('cant do box-cox with negatives')
    logy = log(Yvec)
    xvx = crossprod(Xmat)
    boxcox = optimize(
        logLbc, interval=c(-1.5, 2.5),
        y=Yvec, x=Xmat,
        logy = logy, xvx = xvx, 
        xvxinv=solve(xvx),
        sumlogy = sum(logy)
    )$min
    boxcox = round(boxcox, 1)
    Sboxcox = seq(from= boxcox -0.1*Nboxcox, by=0.1, len=2*Nboxcox+1 )
  } else {
    Sboxcox = boxcox
  }
  
  if(length(Sboxcox) == 1 | Ny > 1) {
    YrepAdd = rep(0, Ny)
  } else {
    YrepAdd = rep(NA, Ny)
    theOnes = abs(Sboxcox-1) < bceps
    logy = log(Yvec)
    sumlogy = sum(logy)
    YrepAdd[theOnes] = 1
    YrepAdd[!theOnes] = -2*(sumlogy)*(
          Sboxcox[!theOnes]-1
          )
    Yorig = Yvec
    Yvec = NULL
    for(D in Sboxcox)
      Yvec = cbind(Yvec, exp(bc*logy)-1)
    
    theZeros = abs(Sboxcox) < bceps
    Yvec[,theZeros] = logy
    YrepAdd[theZeros] = 2*(sumlogy)
    Yvec[,theOnes] = Yorig
  }
  
  

  obsCov = as.matrix(cbind(Yvec, Xmat))

  fromC = .Call('gmrfLik',
      Q, 
      obsCov, 
      1/propNugget, 
      as.double(YrepAdd))

  xyseq = seq(1, length(propNugget)*ncol(obsCov)^2)
  ssq=array(fromC[xyseq], 
      dim=c(ncol(obsCov),ncol(obsCov),length(propNugget)),
      dimnames=list(colnames(obsCov), colnames(obsCov), propNugget))

  ml = array(
      fromC[-xyseq], 
      dim=c(Ny, length(propNugget), 6),
      dimnames=list(
          colnames(obsCov)[1:Ny], 
          propNugget, 
          c('det','detreml','m2logL.ml', 'm2logL.reml', 'varMl', 'varReml'))
  )
  
  pnMat = matrix(
      propNugget, Ny, length(propNugget),
      byrow=TRUE)
  
  ml = abind::abind(ml, tausq.ml=ml[,,'varMl'],
      tausq.reml=ml[,,'varReml'],
      xisq.ml = ml[,,'varMl']/pnMat,
      xisq.reml = ml[,,'varMl']/pnMat,
      oneminuasr = array(oneminusar, dim(ml)[-3]),
      propNugget=pnMat,
      along=3)
  ml[,'0',c('tausq.ml', 'tausq.reml')] = 0
  ml[,'0',c('xisq.ml', 'xisq.reml')] = 
      ml[,'0',c('varMl', 'varReml')]
  
  betaHat = ssq[-(1:Ny),1:Ny,]
  dimnames(betaHat)[[1]] = paste(
      dimnames(betaHat)[[1]], '.betaHat',sep=''
      )
  ml = abind::abind(ml, aperm(betaHat,c(2,3,1)), along=3)    
      
  seBetaHat = apply(ssq[-(1:Ny),-(1:Ny),],3,diag)

  seBetaHat = array(rep(seBetaHat, Ny), 
      dim=c(dim(seBetaHat), Ny),
      dimnames=c(dimnames(seBetaHat), list(colnames(obsCov)[1:Ny])))

  seBetaHat = aperm(seBetaHat, 3:1)
  se.ml = seBetaHat * array(
      rep(ml[,,'varMl'], dim(seBetaHat)[3]),
          dim=c(dim(ml)[-3], dim(seBetaHat)[3])
      )
  dimnames(se.ml)[[3]] = paste(
      dimnames(se.ml)[[3]], '.se.ml', sep=''
      )    
  se.reml = seBetaHat *array(
      rep(ml[,,'varReml'], dim(seBetaHat)[3]),
      dim=c(dim(ml)[-3], dim(seBetaHat)[3])
  )
  dimnames(se.reml)[[3]] = paste(
      dimnames(se.reml)[[3]], '.se.reml', sep=''
  )    
  
  ml = abind::abind(ml, se.ml, se.reml, along=3)

  ml = abind::abind(ml, logL.ml = -ml[,,'m2logL.ml']/2,
      logL.reml = -ml[,,'m2logL.reml']/2,along=3)
  
  res = abind::abind(ml, 
        range = array(thepar$theo['range'],dim=dim(ml)[-3]), 
        shape = array(thepar$theo['shape'],dim=dim(ml)[-3]), 
        shape_optimal = array(thepar$optimal['shape'],dim=dim(ml)[-3]), 
        along=3)
    
    res = aperm(res,c(3,1,2))
    
    res = drop(res)
    

  attributes(res)$Qinfo = thepar
  
  res
  
}
