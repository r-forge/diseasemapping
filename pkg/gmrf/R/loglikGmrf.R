
bceps = 0.0001

logLbc = function(bc, y, x, 
    logy=log(y), 
    xvx = crossprod(x),     
    xvxinv = solve(xvx),
    sumlogy = sum(logy)
){
  
  if(abs(bc)< bceps){
    ybc = logy
  } else if (abs(bc-1)>bceps){
    ybc = (exp(bc*logy)-1)/bc
  } else {
    ybc = y
  }
  twologj = 2*(sumlogy)*(bc-1)
  
  
  betaHat = xvxinv %*% crossprod(x, ybc)
  ssq = ybc - x %*% betaHat
  
  length(y)*log(sum(ssq^2)) - twologj
  
}

loglikGmrfOneRange = function(
    oneminusar,
    Yvec, Xmat, NN, 
    propNugget=0,
    shape=1,
    boxcox=1,
    fixBoxcox=TRUE,
    boxcoxSeq=c(len=11, by=0.01),
    reml=TRUE,
    adjustEdges=FALSE
) {
  
  if(length(propNugget))
    propNugget = sort(unique(c(0, propNugget)))
  
  Q =  geostatsp::maternGmrfPrec(NN,
      param=c(shape=as.vector(shape),
          oneminusar=as.vector(oneminusar[1]),
          conditionalVariance=1),
      adjustEdges=adjustEdges)
  
  thepar=attributes(Q)$param
  
  
  Yvec = as.matrix(Yvec)
  Ny = ncol(Yvec)
  Nobs = nrow(Yvec)
  
  if(!fixBoxcox){
    if(Ny != 1) warning('cant do box-cox with more than one dataset')
    if(any(Yvec<0)) warning('cant do box-cox with negatives')
    logy = log(Yvec)
    xvx = crossprod(Xmat)
    boxcox = optimize(
        logLbc, interval=c(-1.5, 2.5),
        y=Yvec, x=Xmat,
        logy = logy, xvx = xvx, 
        xvxinv=solve(xvx),
        sumlogy = sum(logy)
    )$min
    Nboxcox = ceiling(boxcoxSeq['len']-1)/2
    Sboxcox = seq(from= boxcox - boxcoxSeq['by']*Nboxcox, by=boxcoxSeq['by'], len=2*Nboxcox+1 )
    Sboxcox = round(Sboxcox/boxcoxSeq['by'])*boxcoxSeq['by']
    Sboxcox = sort(unique(Sboxcox))
  } else {
    Sboxcox = boxcox
  }
  
  if(length(Sboxcox) == 1 | Ny > 1) {
    YrepAdd = rep(0, Ny)
    obsCov = as.matrix(cbind(Yvec, Xmat))
    
  } else {
    theOnes = abs(Sboxcox-1) < bceps
    theZeros = abs(Sboxcox) < bceps
    
    logy = log(Yvec)

    Ymat = matrix(NA, nrow(Yvec), length(Sboxcox))
    colnames(Ymat) = as.character(Sboxcox)

    for(D in Sboxcox)
      Ymat[,as.character(D)] = (exp(D*logy)-1)/D
    
    Ymat[,theZeros] = logy
    Ymat[,theOnes] = Yvec
    
    sumlogy = sum(logy)
    
    YrepAdd = 2*(sumlogy)*(
          Sboxcox -1
          )

    Ymat[,theOnes] = Yvec
    obsCov = as.matrix(cbind(Ymat, Xmat))
    Ny = ncol(Ymat)
  }
  
  


  xisqTausq = 1/propNugget 
#  dyn.unload('../src/gmrfLik.so')
#  dyn.load('../src/gmrfLik.so')
  fromC = .Call('gmrfLik',
      Q, 
      obsCov, 
      xisqTausq, 
      as.double(YrepAdd))

  NxysqTausq = length(fromC)/(8*Ny+ncol(obsCov)^2)
      
  xyseq = seq(1, NxysqTausq*ncol(obsCov)^2)
  ssq=array(fromC[xyseq], 
      dim=c(ncol(obsCov),ncol(obsCov),NxysqTausq),
      dimnames=list(colnames(obsCov), colnames(obsCov), NULL))
   
  mlColNames = c('det','detreml','m2logL.ml', 'm2logL.reml', 
      'varMl', 'varReml', 'xisqTausq')
  ml = array(
      fromC[-xyseq], 
      dim=c(Ny, NxysqTausq, length(mlColNames)),
      dimnames=list(
          colnames(obsCov)[1:Ny], 
          NULL, mlColNames
          )
  )
  

  nuggetUsed = !is.na(ml[1,,'det'])
  ml = ml[,nuggetUsed,]
  ssq = ssq[,,nuggetUsed]

  propNugget= 1/ml[1,,'xisqTausq']
  dimnames(ml)[[2]] = dimnames(ssq)[[3]] =
      as.character(propNugget)
  ml = abind::abind(ml, 
      propNugget=matrix(propNugget, nrow=dim(ml)[1], ncol=dim(ml)[2],byrow=TRUE),
      along=3
  )
  
  if(length(Sboxcox)>1){
    

    bestBC = apply(ml[,,'m2logL.ml'], 2, which.min )
    newml = NULL
    newssq = NULL
    
    Scov = c(NA,seq(Ny+1, ncol(obsCov)))
    
    for(Dbc in 1:length(bestBC)){
      
      newml = abind::abind(newml, ml[bestBC[Dbc],Dbc,,drop=FALSE],along=2)
      
      Scov[1] = bestBC[Dbc]
      newssq = abind::abind(newssq, ssq[Scov, Scov,Dbc,drop=FALSE],along=3)
    }
    dimnames(newssq)[[3]] = dimnames(newml)[[2]] = dimnames(ml)[[2]]
    dimnames(newml)[[1]] = 'optbc'
    newml = abind::abind(newml, boxcox=matrix(Sboxcox[bestBC], nrow=1), along=3)
    ssq=newssq
    ml=newml
  }
  

  newml = abind::abind(ml, 
      tausq.ml=ml[,,'varMl',drop=FALSE],
      tausq.reml=ml[,,'varReml',drop=FALSE],
      xisq.ml = ml[,,'varMl',drop=FALSE]/ml[,,'propNugget',drop=FALSE],
      xisq.reml = ml[,,'varMl',drop=FALSE]/ml[,,'propNugget',drop=FALSE],
      oneminuasr = array(oneminusar, dim(ml)[-3]),
      along=3)
  dimnames(newml)[[3]][-seq(1,dim(ml)[3])] = 
      c('tausq.ml','tausq.reml','xisq.ml','xisq.reml','oneminusar')
  ml = newml
  ml[,'0',c('tausq.ml', 'tausq.reml')] = 0
  ml[,'0',c('xisq.ml', 'xisq.reml')] = 
      ml[,'0',c('varMl', 'varReml')]
  
  Ny = dim(ml)[1]
  betaHat = ssq[-(1:Ny),1:Ny,,drop=FALSE]
  dimnames(betaHat)[[1]] = paste(
      dimnames(betaHat)[[1]], '.betaHat',sep=''
      )
  ml = abind::abind(ml, aperm(betaHat,c(2,3,1)), along=3)    
      
  seBetaHat = apply(ssq[-(1:Ny),-(1:Ny),,drop=FALSE],3,diag)

  seBetaHat = array(rep(seBetaHat, Ny), 
      dim=c(dim(seBetaHat), Ny),
      dimnames=c(dimnames(seBetaHat), list(colnames(obsCov)[1:Ny])))

  seBetaHat = aperm(seBetaHat, 3:1)
  se.ml = seBetaHat * array(
      rep(ml[,,'varMl',drop=FALSE], dim(seBetaHat)[3]),
          dim=c(dim(ml)[-3], dim(seBetaHat)[3])
      )
  dimnames(se.ml)[[3]] = paste(
      dimnames(se.ml)[[3]], '.se.ml', sep=''
      )    
  se.reml = seBetaHat *array(
      rep(ml[,,'varReml',drop=FALSE], dim(seBetaHat)[3]),
      dim=c(dim(ml)[-3], dim(seBetaHat)[3])
  )
  dimnames(se.reml)[[3]] = paste(
      dimnames(se.reml)[[3]], '.se.reml', sep=''
  )    
  
  ml = abind::abind(ml, se.ml, se.reml, along=3)
  
  logL.ml = -ml[,,'m2logL.ml',drop=FALSE]/2
  logL.reml = -ml[,,'m2logL.reml',drop=FALSE]/2
  dimnames(logL.reml) = dimnames(logL.ml) = NULL
  ml = abind::abind(ml, 
      logL.ml = logL.ml,
      logL.reml = logL.reml ,along=3)
  
  
  res = abind::abind(ml, 
        range = array(thepar$theo['range'],dim=dim(ml)[-3]), 
        shape = array(thepar$theo['shape'],dim=dim(ml)[-3]), 
        shape_optimal = array(thepar$optimal['shape'],dim=dim(ml)[-3]), 
        along=3)
    
    res = aperm(res,c(3,1,2))
    
    res = drop(res)
    if(is.matrix(res)){
      res = t(res)
      if(reml){
        theMin = which.min(res[,'m2logL.reml'])
      } else{
        theMin = which.min(res[,'m2logL.ml'])
        
      }
      res = rbind(opt=res[theMin,], res)
    }

  attributes(res)$Qinfo = thepar
  
  res
  
}
