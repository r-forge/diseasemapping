#+ setup
library('geostatsp')
#'

#' # simulated data


#+ simData
myRaster = squareRaster(extent(0,8000,0,6000), 50)
myParam=c(oneminusar=0.1, conditionalVariance=100,shape=2)
myQ = maternGmrfPrec(myRaster, param=myParam)
attributes(myQ)$info$optimalShape
set.seed(0)
mySim = RFsimulate(attributes(myQ)$info$optimalShape, myRaster)

otherPar = c(intercept=1, beta = 2, tau=10)
myCov = myRaster
values(myCov) = rep(seq(-1,1,len=ncol(myCov)), nrow(myCov))

myLambda = 1 + 2 * myCov + mySim
myY = myLambda
values(myY)=rnorm(prod(dim(myLambda)),values(myLambda), sd=10) 

names(myCov) = 'x'
names(myY) = gsub("^layer\\.","sim", names(mySim))
#'





#' grid search	

#+ simGrid
myResR = lgm(formula = sim ~ x, 
  data=raster::stack(myY, myCov), 
  oneminusar = seq(0.2, 0.7,len=25),
  nugget = seq(0, 0.001,len=21), shape=2, 
  adjustEdges=TRUE,
  mc.cores=1+(.Platform$OS.type=='unix') )		
#'

#+ simPlot
Sbreaks = c(-100000,-50,-20,-10, -5,  -2, -1,0)

myCol = mapmisc::colourScale(
  breaks = Sbreaks + max(myResR$array[,-1,'logLreml',]),
  style='fixed',
  col=terrain.colors
)

image(	
  myResR$array[1,-1,'propNugget',1], 
  myResR$array[1,1,'oneminusar',], 
  myResR$array[1,-1,'logLreml',],
  xlab = 'propNugget', ylab='oneminusar',
col=myCol$col, breaks=myCol$breaks)
mapmisc::legendBreaks("topright", breaks = Sbreaks, col=myCol$col)


points(myResR$param['propNugget'], myResR$param['oneminusar'])
#'


#+ checkResults
myResR = lgm(formula = sim ~ x, 
  data=raster::stack(myY, myCov), 
  oneminusar = seq(0.2, 0.7,len=25),
  nugget = seq(0, 0.001,len=21), shape=2, 
  adjustEdges=TRUE,
  mc.cores=1+(.Platform$OS.type=='unix') )		
oneRes =myResR$array[,3,,3]

Q = maternGmrfPrec(
  myRaster, 
  param=c(conditionalVariance=1, oneRes[c('oneminusar','shape')]),
  adjustEdges=TRUE)

Qinv = as.matrix(solve(Q))
V =  oneRes['propNugget'] * Qinv + Diagonal(nrow(Q))
Vinv = solve(V)

obsCov = cbind(y=values(myY), intercept=1, x=values(myCov))

(xProdQ =  crossprod(obsCov, Vinv) %*% obsCov)

(betahat= drop(solve(xProdQ[Xseq,Xseq]) %*% xProdQ[Xseq, -Xseq]))

myResR$array[,
  as.character(1/oneRes['propNugget']), 
  c('(Intercept)BetaHat','xBetaHat'),
  as.character(oneRes['oneminusar'])]

resid = obsCov[,1] - betahat[1] - obsCov[,3]*betahat[2]
(varhat = diag(crossprod(resid, Vinv) %*% resid)/Nobs)
myResR$array[,
  as.character(1/oneRes['propNugget']), 
  'tausqHatMl',
  as.character(oneRes['oneminusar'])]




xisqTausq = c(0, myResR$array[,-1,'propNugget',1])
myResR$array[,,'propNugget',1]
jacobian=0



Nxy = ncol(obsCov)
Ny = length(jacobian)
Nobs = nrow(obsCov)
NxisqTausq = length(xisqTausq)
Nxysq = Nxy^2
logLstart = NxisqTausq*Nxysq*2
Lseq = 1:logLstart
mlColNames = c('det','detReml','m2logLml', 
  'm2logLreml', 'profiledVarianceHatMl',
  'profiledVarianceHatReml', 'xisqTausq')

mlDim = c(y=Ny,
  varRatio=NxisqTausq, 
  output=length(mlColNames))

fromC = .Call('gmrfLik',
  Q, 
  obsCov, 
  as.double(xisqTausq), 
  reml=TRUE,
  as.double(jacobian),
  as.double(c(-999,0,0))
)

stuff3 = loglikGmrf(
  Yvec=obsCov[,1], Xmat=obsCov[,-1], 
  NN = NNmat(myRaster), 
  oneminusar = oneRes['oneminusar'],
  propNugget= 1/xisqTausq[-1],
  boxcox=1,
  fixBoxcox=TRUE,
  seqBoxcox=outer(c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5), c(-1,1)),
  shape=oneRes['shape'],
  reml=TRUE,
  adjustEdges=TRUE,
  jacobian = 0,
  control=list(
    oneminusar = c(lower=1e-3, upper=0.8, tol= .Machine$double.eps^0.25),
    propNugget = c(lower=1e-3, upper=1e3, tol= .Machine$double.eps^0.25)
  ),
  mc.cores=4
) 

Xseq = 2:3

stuff = array(
  as.vector(fromC[
      seq(logLstart+1, len=prod(mlDim))
    ]), 
  dim=mlDim,
  dimnames=list(
    colnames(obsCov)[1:Ny], 
    as.character(xisqTausq), 
    mlColNames
  ))

stuff2 = array(fromC[Lseq], 
  dim=c(Nxy, Nxy, NxisqTausq,2),
  dimnames = list(
    colnames(obsCov), colnames(obsCov), 
    as.character(xisqTausq),
    c('ssq','beta')
  ))

stuff[1,1,'det']
determinant(Q)$modulus

stuff[1,1,'detReml']
determinant((crossprod(obsCov, Q)%*% obsCov)[Xseq,Xseq])$modulus

stuff2[,,1,1]
(xProdQ =  crossprod(obsCov, Q) %*% obsCov)

stuff2[Xseq,1,1,2]
(betahat=drop(solve(xProdQ[Xseq,Xseq]) %*% xProdQ[Xseq, -Xseq]))

resid = obsCov[,1] - betahat[1] - obsCov[,3]*betahat[2]
(varhat = diag(crossprod(resid, Q) %*% resid)/Nobs)
stuff[,1,'profiledVarianceHatMl']
stuff3$extras$ml[,'propNugget=0','profiledVarianceHatMl',as.character(oneRes['oneminusar'])]
myResR$model$extras$ml[,'propNugget=0','profiledVarianceHatMl',as.character(oneRes['oneminusar'])]
myResR$array[,
  '0', 
  'tausqHatMl',
  as.character(oneRes['oneminusar'])]


stuff[,,'m2logLml']
Nobs*log(2*pi) + Nobs*log(varhat) - determinant(Q)$modulus + stuff2[-Xseq,-Xseq,1,'beta']/varhat
Nobs + Nobs*log(2*pi) - determinant(Q)$modulus + Nobs*log(varhat) 
myResR$model$extras$ml[,'propNugget=0','m2logLml', as.character(oneRes['oneminusar'])]

-2*mvtnorm::dmvnorm(
  resid, 
  sigma = stuff[,1,'profiledVarianceHatMl'] * Qinv, log=TRUE)

betahat
myResR$array[,'0',c('logLml','mlDeviance','(Intercept)BetaHat','xBetaHat'), as.character(oneRes['oneminusar'])]

# with nugget

stuff2[,,2,1]
(xProdQ =  crossprod(obsCov, Vinv) %*% obsCov)

stuff2[Xseq,1,2,2]
(betahat= drop(solve(xProdQ[Xseq,Xseq]) %*% xProdQ[Xseq, -Xseq]))
stuff3$extras$ml[,
  paste('propNugget=', oneRes['propNugget'],sep=''),
  'profiledVarianceHatMl',as.character(oneRes['oneminusar'])]

myResR$array[,as.character(1/oneRes['propNugget']),, as.character(oneRes['oneminusar'])]


resid = obsCov[,1] - betahat[1] - obsCov[,3]*betahat[2]
(varhat = diag(crossprod(resid, Vinv) %*% resid)/Nobs)
stuff[,2,'profiledVarianceHatMl']
myResR$model$extras$ml[,
  paste('propNugget=',oneRes['propNugget'],sep=''),,
  as.character(oneRes['oneminusar'])]
myResR$model$extras$ml[,
  3,,
  as.character(oneRes['oneminusar'])]

stuff[,,'m2logLml']
Nobs*log(2*pi) + Nobs*log(varhat) - determinant(Vinv)$modulus + stuff2[-Xseq,-Xseq,2,'beta']/varhat
Nobs + Nobs*log(2*pi) - determinant(Vinv)$modulus + Nobs*log(varhat)

V =  oneRes['propNugget'] * Qinv + Diagonal(nrow(Q))

-2*mvtnorm::dmvnorm(resid, sigma = varhat*as.matrix(V), log=TRUE)




myResR$model$extras$ml[,1:2,c('det','detReml'),   as.character(oneRes['oneminusar'])]





Vinv = solve(V)


myResR$model$extras$ssq[,,
  as.character(1/oneRes['propNugget']),
  ,
  as.character(oneRes['oneminusar'])]

crossprod(obsCov, Vinv) %*% obsCov


crossprod(obsCov) - SxisqtausqFull[2] * 
  t(obsCov) %*% QalmostInv %*% obsCov


varY = as.matrix(oneRes['xisqHatMl']*solve(qTest) + oneRes['tausqHatMl']*diag(ncol(qTest)))

myResR$model$extras$ml[,
  paste('propNugget=',oneRes['propNugget'],sep=''),,
  as.character(oneRes['oneminusar'])]

determinant(varY)$modulus

Yresid = myY - oneRes['(Intercept)BetaHat'] - oneRes['xBetaHat'] * myCov

mvtnorm::dmvnorm(values(Yresid), sigma = varY, log=TRUE)

oneRes['logLml']
#' 

# swiss rain

data('swissRainR')

anotherx = raster(swissRainR[['alt']])
values(anotherx) = seq(0,1,len=ncell(anotherx))
names(anotherx) = "myvar"

swissRainR2 = brick(swissRainR[['alt']], 
  sqrt(swissRainR[['prec1']]),
  anotherx)


swissResR =  lgm(
  formula=layer ~ alt+ myvar, 
  data=swissRainR2, shape=2,
  oneminusar = exp(seq(log(0.001), log(0.025), len=11)),
  nugget = exp(seq(log(100), log(50000), len=11)),
  adjustEdges=TRUE,
  mc.cores=1+(.Platform$OS.type=='unix') )		


myCol = mapmisc::colourScale(
    breaks = Sbreaks + max(swissResR$array[,-1,'logLreml',]),
    style='fixed',
    col=terrain.colors
)

image(	
    swissResR$array[1,-1,'propNugget',1], 
    swissResR$array[1,1,'oneminusar',], 
    swissResR$array[1,-1,'logLreml',],
    xlab = 'propNugget', ylab='oneminusar',
    log='xy',
    col=myCol$col, breaks=myCol$breaks)
mapmisc::legendBreaks("topright", breaks = Sbreaks, col=myCol$col)



if(Sys.info()['user'] =='patrick' & FALSE) {
  
  
  
# boxcox
  yBC = sqrt(myY + 1 - minValue(myY))
  names(yBC) = names(myY)
  myResBC = lgm(
    formula = sim ~ x, 
    data=raster::stack(yBC, myCov), 
    oneminusar = seq(0.02, 0.3,len=24),
    nugget = seq(0, 2,len=40), 
    shape=2, 
    mc.cores=1+(.Platform$OS.type=='unix'), 
    fixBoxcox=FALSE,
    adjustEdges=FALSE)
  
  
  if(!interactive()) pdf("profLboxcox.pdf")
  plot(myResBC$profL$boxcox,type='o', ylim=max(myResBC$profL$boxcox[,2])-c(3,0))
  if(!interactive()) dev.off()
  
  myResBC$param
  
  myCol = mapmisc::colourScale(
    breaks = Sbreaks,
    style='fixed',
    col=terrain.colors
  )
  
  if(!interactive()) pdf("profLwithboxcox.pdf")
  image(	myResBC$array[,-1,'propNugget',1], 
    myResBC$array[,1,'oneminusar',], 
    myResBC$array[,-1,'logLreml',],
    col=myCol$col, breaks=myCol$breaks+max(myResBC$array[,,'logLreml',]))
  mapmisc::legendBreaks("topright",  myCol)
  points(myResBC$param['propNugget'], myResBC$param['oneminusar'])
  if(!interactive()) dev.off()
  
# optimizing, doesn't work
  
# optimize propNugget
  myResRopt = lgm(
    formula = sim ~ x, 
    data=raster::stack(myY, myCov), 
    oneminusar = seq(0.05, 0.2,len=12),
    shape=2)		

  if(!interactive()) pdf("doesntwork.pdf")
  plot(myResRopt$array[,,'oneminusar',], myResRopt$array[,,'propNugget',])
  
  if(!interactive()) dev.off()	

swissResRoptAr =  lgm(
  formula=layer ~ alt+ myvar, 
  data=swissRainR2, shape=2,
  oneminusar = seq(0.1, 0.5, len=6),
  adjustEdges=FALSE
)

swissResRopt =  lgm(
  formula=layer ~ alt+ myvar, 
  data=swissRainR2, shape=2,
  adjustEdges=FALSE
)


swissResRopt$summary

# with edge correction.  
# time consuming, only run this if Patrick is checking



# optimize only nugget
swissResROptNug =  lgm(
  formula=layer ~ alt+ myvar, 
  data=swissRainR2, shape=2,
  oneminusar=seq(0.05, 0.1, len=12),
  adjustEdges=FALSE,fixNugget=TRUE,
   mc.cores=1+(.Platform$OS.type=='unix')
)

plot(swissResROptNug$profL$range, type='l')
}



