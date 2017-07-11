#+ setup
library('geostatsp')
#'

#' # simulated data


#+ simData
myRaster = squareRaster(extent(0,6000,0,6000), 60)
myParam=c(oneminusar=0.1, conditionalVariance=2.5^2,shape=2)
myQ = maternGmrfPrec(myRaster, param=myParam)
attributes(myQ)$info$optimalShape
set.seed(0)

mySimGRF = RFsimulate(attributes(myQ)$info$optimalShape, myRaster)
mySimGmrf = raster(myRaster)
cholQ = chol(myQ)#, perm=FALSE, ldl=FALSE)
set.seed(0)
values(mySimGmrf) = as.vector(solve(t(cholQ), rnorm(ncell(myRaster))))
plot(mySimGmrf)


values(mySimGmrf) = as.vector(Matrix::solve(cholQ, solve1, system='P'))

otherPar = c(intercept=1, beta = 2, tau=10)
myCov = myRaster
values(myCov) = rep(seq(-1,1,len=ncol(myCov)), nrow(myCov))

myLambda = otherPar['intercept'] + otherPar['beta'] * myCov + mySim
myY = myLambda
values(myY)=rnorm(prod(dim(myLambda)),values(myLambda), sd=otherPar['tau']) 

names(myCov) = 'x'
names(myY) = gsub("^layer\\.","sim", names(mySim))
#'



#+ plotSim, fig.cap='simulated data', fig.subcap =c('random effect','observed')
plot(mySim)

plot(myY)
#'

#' grid search	

#+ simGrid
myResR = lgm(formula = sim ~ x, 
    data=raster::stack(myY, myCov), 
    oneminusar = exp(seq(log(0.01), log(0.1),len=25)),
    nugget = exp(seq(log(10), log(1000),len=21)), shape=2, 
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
    as.numeric(dimnames(myResR$array)$propNugget)[-1], 
     as.numeric(dimnames(myResR$array)$oneminusar),
     myResR$array[1,-1,'logLreml',],
    xlab = 'propNugget', ylab='oneminusar',
    log='xy',
    col=myCol$col, breaks=myCol$breaks)
mapmisc::legendBreaks("topright", breaks = Sbreaks, col=myCol$col)


points(myResR$param['propNugget'], myResR$param['oneminusar'])
#'


#+ results
myResR$param
attributes(myQ)$info$optimalShape
otherPar

#'

#' checkResults

#+ checkResultsSetup

oneRes = c(
    oneminusar = as.numeric(dimnames(myResR$array)$oneminusar[5]), 
    shape=2, 
    propNugget = as.numeric(dimnames(myResR$array)$propNugget[4]))

Q = maternGmrfPrec(
    myRaster, 
    param=c(conditionalVariance=1, oneRes[c('oneminusar','shape')]),
    adjustEdges=TRUE)

Qinv = as.matrix(solve(Q))

obsCov = cbind(y=values(myY), intercept=1, x=values(myCov))
Xseq = 2:ncol(obsCov)
Nobs = nrow(obsCov)
#'

#' no nugget

#+ checkNoNugget
(xProdQ =  crossprod(obsCov, Q) %*% obsCov)
xProdQ - myResR$model$extras$ssq[,,1,'ssq',  as.character(oneRes['oneminusar'])]

myResR$array[,1,c('(Intercept)BetaHat','xBetaHat'),as.character(oneRes['oneminusar'])]
myResR$model$extras$ssq[Xseq,1,1,'beta',  as.character(oneRes['oneminusar'])]
(betahat=drop(solve(xProdQ[Xseq,Xseq]) %*% xProdQ[Xseq, -Xseq]))

myResR$array[,1,c('xisqHatMl','xisqHatReml'),as.character(oneRes['oneminusar'])]
myResR$model$extras$ml[,1,c('profiledVarianceHatMl','profiledVarianceHatReml'),  as.character(oneRes['oneminusar'])]
resid = obsCov[,1] - betahat[1] - obsCov[,3]*betahat[2]
(varhat = diag(crossprod(resid, Q) %*% resid)/(nrow(obsCov) - c(ml=0,reml=length(Xseq))))


c(
    -0.5*myResR$model$extras$ml[,'0','m2logLml', as.character(oneRes['oneminusar'])],
    myResR$array[,1,'logLml',as.character(oneRes['oneminusar'])],
    mvtnorm::dmvnorm(
        resid, 
        sigma = varhat['ml'] * Qinv, log=TRUE)
)
#'

#' with nugget

#+ checkWithNugget
V =  (1/oneRes['propNugget']) * Qinv + Diagonal(nrow(Q))
Vinv = solve(V)

(xProdQ =  crossprod(obsCov, Vinv) %*% obsCov)
myResR$model$extras$ssq[,,  as.character(oneRes['propNugget']),'ssq',  as.character(oneRes['oneminusar'])]


(betahat= drop(solve(xProdQ[Xseq,Xseq]) %*% xProdQ[Xseq, -Xseq]))
myResR$array[,
    as.character(oneRes['propNugget']), 
    c('(Intercept)BetaHat','xBetaHat'),
    as.character(oneRes['oneminusar'])]

resid = obsCov[,1] - betahat[1] - obsCov[,3]*betahat[2]
(varhat = diag(crossprod(resid, Vinv) %*% resid)/(nrow(obsCov) - c(ml=0,reml=length(Xseq))))
myResR$array[,
    as.character(oneRes['propNugget']), 
    c('tausqHatMl','tausqHatReml'),
    as.character(oneRes['oneminusar'])]

c(
    -0.5*myResR$model$extras$ml[,as.character(oneRes['propNugget']),'m2logLml', as.character(oneRes['oneminusar'])],
    myResR$array[,as.character(oneRes['propNugget']),'logLml',as.character(oneRes['oneminusar'])],
    mvtnorm::dmvnorm(
        resid, 
        sigma = as.matrix(varhat['ml'] * V), log=TRUE)
)
#' 

#' # swiss rain

#+ swissRain
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
#'



#' boxcox
#+ boxCox
if(Sys.info()['user'] =='patrick' & FALSE) {
  
  
  yBC = sqrt(myY + 1 - minValue(myY))
  names(yBC) = names(myY)
  myResBC = lgm(
      formula = sim ~ x, 
      data=raster::stack(yBC, myCov), 
      oneminusar = seq(0.1,0.5,len=24),
      nugget = seq(1, 10,len=40), 
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
  image(	
      myResBC$profL$twoDim$x, 
      myResBC$profL$twoDim$y,
      myResBC$profL$twoDim$z,
      ylab = 'range', xlab='propNugget',
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



#'