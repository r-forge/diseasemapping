library('geostatsp')


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



# simulated data


# grid search	


myResR = lgm(formula = sim ~ x, 
  data=raster::stack(myY, myCov), 
  oneminusar = seq(0.2, 0.7,len=25),
  nugget = seq(0, 0.001,len=21), shape=2, 
  adjustEdges=TRUE,
  mc.cores=1+(.Platform$OS.type=='unix') )		

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



