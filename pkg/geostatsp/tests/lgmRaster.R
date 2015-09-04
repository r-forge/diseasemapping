library('geostatsp')


myRaster = squareRaster(extent(0,8000,0,6000), 50)
myParam=c(oneminusar=0.1, conditionalVariance=100,shape=2)
myQ = maternGmrfPrec(myRaster, param=myParam)

set.seed(0)
mySim = RFsimulate(attributes(myQ)$param$theo, myRaster)

otherPar = c(intercept=1, beta = 2, tau=10)
myCov = myRaster
values(myCov) = rep(seq(-1,1,len=ncol(myCov)), nrow(myCov))



myLambda = 1 + 2 * myCov + mySim
myY = myLambda
values(myY)=rnorm(prod(dim(myLambda)),values(myLambda), sd=10) 

names(myCov) = 'x'
names(myY) = gsub("^layer\\.","sim", names(mySim))

if(Sys.info()['user'] =='patrick') {

# grid search	
myResR = lgm(formula = sim ~ x, 
  	data=raster::stack(myY, myCov), 
  	oneminusar = seq(0.05, 0.2,len=12),
  	nugget = seq(0.25, 4,len=12), shape=2, mc.cores=1)		

image(	myResR$array[,-1,'propNugget',1], 
		myResR$array[,1,'oneminusar',], 
		myResR$array[,-1,'logLreml',])

# boxcox
	yBC = sqrt(myY + 50)
	names(yBC) = names(myY)
	myResBC = lgm(formula = sim ~ x, 
  		data=raster::stack(yBC, myCov), 
  		oneminusar = seq(0.05, 0.2,len=12),
  		nugget = seq(0.25, 4,len=12), shape=2, mc.cores=1, 
			fixBoxcox=FALSE,
			adjustEdges=FALSE)
	
	image(	myResR$array[,-1,'propNugget',1], 
			myResR$array[,1,'oneminusar',], 
			myResR$array[,-1,'logLreml',])
	

# optimize propNugget
	myResRopt = lgm(
			formula = sim ~ x, 
  		data=raster::stack(myY, myCov), 
  		oneminusar = seq(0.05, 0.2,len=12),
			shape=2)		
	pdf("doesntwork.pdf")
	plot(myResRopt$array[,,'oneminusar',], myResRopt$array[,,'propNugget',])
	dev.off()	
}






data("swissRainR")

anotherx = raster(swissRainR[['alt']])
values(anotherx) = seq(0,1,len=ncell(anotherx))
names(anotherx) = "myvar"

swissRainR2 = brick(swissRainR[['alt']], 
		sqrt(swissRainR[['prec1']]),
		anotherx)

if(Sys.info()['user'] =='patrick') {
	
swissResR =  lgm(
    formula=layer ~ alt+ myvar, 
    data=swissRainR2, shape=2,
		oneminusar = seq(0.01, 0.1, len=12),
		nugget = seq(0, 1, len=20),
    adjustEdges=TRUE
)

image(	swissResR$array[,,'propNugget',1], 
		swissResR$array[,1,'oneminusar',], 
		swissResR$array[,,'logLreml',])
}

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
if(Sys.info()['user'] =='patrick') {


# optimize only nugget
swissResROptNug =  lgm(
    formula=layer ~ alt+ myvar, 
    data=swissRainR2, shape=2,
    oneminusar=seq(0.05, 0.1, len=12),
    adjustEdges=FALSE,fixNugget=TRUE,
    mc.cores=c(1,2)[1+(.Platform$OS.type=='unix')]
)

plot(swissResROptNug$profL$range, type='l')

}


