library(geostatsp)

myraster = raster(extent(0,8000,0,6000), ncols=40,nrows=30)
myrasterBig = extend(myraster, 
		extend(extent(myraster),6*xres(myraster))
)

theNN = NNmat(myraster)
themodel = c(range=5*xres(myraster),shape=2,variance=900)
if(T){
	theU = RFsimulate(myrasterBig,model=themodel)
	theU = raster::crop(theU, myraster)
} else {
	thePrec = maternGmrfPrec(theNN, 
		param=c(themodel,cellSize=xres(myraster)),
		adjustEdges=TRUE,adjustParam=TRUE)

thePrec2 = maternGmrfPrec(theNN, 
		param=c(themodel,cellSize=xres(myraster)),
		adjustEdges=TRUE,adjustParam=FALSE)

	theVar = matern(myraster,param=themodel)
	theprod = thePrec %*% theVar
	par(mfrow=c(2,2))
	hist(diag(theprod),breaks=60)
	hist(theprod[lower.tri(theprod,diag=F)],breaks=60)
	
	theCholP = chol(thePrec)
	theCholP2 = chol(thePrec2)
	theCholV = chol(theVar)
	
	theprod = theCholV %*% theCholP2
	hist(diag(theprod),breaks=60)
	hist(theprod[upper.tri(theprod,diag=F)],breaks=60)
	
	
	resM = NULL
	myfun = function(qq) c(m=mean(qq), v=var(qq))
	for(D in 1:40) {
	theZ = rnorm(ncell(myraster))
	resM = cbind(c(
					p=myfun(as.vector(solve(theCholP,theZ))),
					v=myfun(as.vector(theCholV%*%theZ)),
					p2=myfun(as.vector(solve(theCholP2,theZ)))
					),resM)	
}

apply(resM,1, quantile)
apply(resM,1, sd)

	theU =theUP= myraster
	values(theUP) = as.vector(solve(theCholP, theZ))
	values(theU) = as.vector(theCholV %*% theZ)
	par(mfrow=c(2,2))
	plot(theU)
	plot(theUP)
	plot(theU/theUP)
	plot(theU-theUP)
}


if(FALSE) {
	
	myrasterH = raster(extent(0,40*1000,0,30*1000), 
			ncols=ncol(myraster),nrows=nrow(myraster))
	res(myrasterH)
	theNNH = NNmat(myrasterH)
	
	paramH=c(range=4*xres(myrasterH),
			shape=2,var=100,cellSize=xres(myrasterH))
	paramH2 = c(range=4,
			shape=2,var=100)
	
	varMat = matern(myrasterH,param=paramH)
	
	precNeither = maternGmrfPrec(theNNH,param=paramH2,
			adjustEdges=FALSE,adjustParam=FALSE)
	precEdge = maternGmrfPrec(theNNH,param=paramH2,
			adjustEdges=TRUE,adjustParam=FALSE)

	precAdj = maternGmrfPrec(theNNH,param=paramH2,
			adjustEdges=FALSE,adjustParam=TRUE,adjustShape=FALSE)
	
	precBoth = maternGmrfPrec(theNNH,param=paramH2,
			adjustEdges=TRUE,adjustParam=TRUE,adjustShape=FALSE)
	
	theZ = rnorm(ncell(myrasterH))
	resM = NULL
	for(D in c('Neither','Edge','Adj','Both')){
		theChol = chol(get(paste('prec',D,sep='')))
		resM = c(resM,
				mean(as.vector(solve(theChol, theZ))))		
	}
	resM
	
	prodNeither = precNeither %*% varMat
	prodEdge = precEdge %*% varMat
	prodAdj = precAdj %*% varMat
	prodBoth = precBoth %*% varMat

	
	par(mfrow=c(4,2))
	
	hist(Matrix::diag(prodNeither),breaks=60,ylab='neither')
	abline(v=1)
	hist(prodNeither[lower.tri(prodNeither,diag=FALSE)],breaks=60)	
	abline(v=0)
	hist(Matrix::diag(prodEdge),breaks=60,ylab='edge')
	abline(v=1)
	hist(prodEdge[lower.tri(prodEdge,diag=FALSE)],breaks=60)	
	abline(v=0)
	
	
	hist(Matrix::diag(prodAdj),breaks=60,ylab='adj')
	abline(v=1)
	hist(prodAdj[lower.tri(prodAdj,diag=FALSE)],breaks=60)	
	abline(v=0)
	

	hist(Matrix::diag(prodBoth),breaks=60,ylab='both')
	abline(v=1)
	hist(prodBoth[lower.tri(prodBoth,diag=FALSE)],breaks=60)	
	abline(v=0)

	
}


thecov = myraster
values(thecov) = c(rep(1,ncell(thecov)/2),
		rep(4,ncell(thecov)/2))
beta.x=5
theY = theU + beta.x*thecov

theData = brick(theY,thecov)
names(theData) = c("y","x")

thedf = as.data.frame(theData)

Xmat=cbind(inter=1,x=thedf$x)
maternShape=as.numeric(themodel['shape'])

Sar = exp(seq(log(0.025),log(0.25),len=24))


#oneminusar=0.2;shape=maternShape;Yvec=thedf$y;

source("../R/loglikGmrf.R")

resNoNuggetVanilla = loglikGmrf(
		oneminusar=Sar,
		Yvec=thedf$y,
		Xmat=Xmat,
		shape=maternShape,
		NN=theNN,mc.cores=4
		)

resNoNuggetEdge = loglikGmrf(
		oneminusar=Sar,
				Yvec=thedf$y,
				Xmat=Xmat,
				shape=maternShape,
				NN=theNN,mc.cores=4,
				adjustEdges=TRUE,	adjustParam=FALSE
		)
#		
resNoNuggetAdj = loglikGmrf(
		oneminusar=Sar,
				Yvec=thedf$y,
				Xmat=Xmat,
				shape=maternShape,
				NN=theNN,mc.cores=4,
				adjustParam=TRUE,adjustEdges=FALSE,
				adjustShape=FALSE
		)
		

resNoNuggetBoth = loglikGmrf(
		oneminusar=Sar,
				Yvec=thedf$y,
				Xmat=Xmat,
				shape=maternShape,
				NN=theNN,mc.cores=4,
				adjustEdges=TRUE,	adjustParam=TRUE,
				adjustShape=FALSE
		)

resNoNuggetShape = loglikGmrf(
			oneminusar=Sar,
				Yvec=thedf$y,
				Xmat=Xmat,
				shape=maternShape,
				NN=theNN,mc.cores=4,
				adjustEdges=TRUE,	adjustParam=TRUE,
				adjustShape=TRUE
		)
		
		
		col = c(vanilla = 'black', adj='green',edge='blue',both='red',shape='orange')
		
		
		scale = sqrt(8*themodel['shape'])*xres(myraster)/
				themodel["range"]
		a = (scale^2 + 4) 
		oneminusar = 1-4/a
		
		par(mfrow=c(2,2))		
		D ='logL.ml' 
		
		plot(resNoNuggetEdge['oneminusar',],
				resNoNuggetEdge[D,],
				xlab='1-ar',ylab=D,
	 			ylim=c(median(resNoNuggetVanilla[D,]),
	 					max(resNoNuggetShape[D,])),
				col=col['edge'],type='o',log='x',lwd=1)
		legend('left',fill=col,legend=names(col))
		lines(resNoNuggetVanilla['oneminusar',],
				resNoNuggetVanilla[D,],
				type='o', col=col['vanilla'],lwd=3)
		lines(resNoNuggetAdj['oneminusar',],
				resNoNuggetAdj[D,],
				type='o', col=col['adj'])		
		lines(resNoNuggetBoth['oneminusar',],
				resNoNuggetBoth[D,],
				type='o', col=col['both'])		
		lines(resNoNuggetShape['oneminusar',],
				resNoNuggetShape[D,],
				type='o', col=col['shape'])		
		

		abline(v= oneminusar)

plot(resNoNuggetEdge['rangeInCells',]*xres(myraster),
		resNoNuggetEdge[D,],
		ylim=c(median(resNoNuggetVanilla[D,]),
				max(resNoNuggetShape[D,])),
				xlab='range',ylab=D,
				col=col['edge'],type='o',log='x',
				lwd=1)

lines(resNoNuggetVanilla['rangeInCells',]*xres(myraster),
		resNoNuggetVanilla[D,],
		type='o', col=col['vanilla'],lwd=3)

lines(resNoNuggetAdj['rangeInCells',]*xres(myraster),
		resNoNuggetAdj[D,],
		type='o', col=col['adj'])		
lines(resNoNuggetBoth['rangeInCells',]*xres(myraster),
		resNoNuggetBoth[D,],
		type='o', col=col['both'])		
lines(resNoNuggetShape['rangeInCells',]*xres(myraster),
		resNoNuggetShape[D,],
		type='o', col=col['shape'])		
abline(v=themodel['range'])
 

Dx='sigmasq.ml'



plot(resNoNuggetEdge[Dx,],
		resNoNuggetEdge[D,],
		ylim=c(median(resNoNuggetVanilla[D,]),
				max(resNoNuggetShape[D,])),
		xlab=Dx,ylab=D,
		col=col['edge'],type='o',log='x',
		lwd=1)

lines(resNoNuggetVanilla[Dx,],
		resNoNuggetVanilla[D,],
		type='o', col=col['vanilla'],lwd=3)

lines(resNoNuggetAdj[Dx,],
		resNoNuggetAdj[D,],
		type='o', col=col['adj'])		
lines(resNoNuggetBoth[Dx,],
		resNoNuggetBoth[D,],
		type='o', col=col['both'])		
lines(resNoNuggetShape[Dx,],
		resNoNuggetShape[D,],
		type='o', col=col['shape'])		
abline(v=themodel['variance'])


Dx='rangeInCells'
D ='logL.reml' 


plot(resNoNuggetEdge[Dx,],
		resNoNuggetEdge[D,],
		ylim=c(median(resNoNuggetVanilla[D,]),
				max(resNoNuggetShape[D,])),
		xlab=Dx,ylab=D,
		col=col['edge'],type='o',log='x',
		lwd=1)

lines(resNoNuggetVanilla[Dx,],
		resNoNuggetVanilla[D,],
		type='o', col=col['vanilla'],lwd=3)

lines(resNoNuggetAdj[Dx,],
		resNoNuggetAdj[D,],
		type='o', col=col['adj'])		
lines(resNoNuggetBoth[Dx,],
		resNoNuggetBoth[D,],
		type='o', col=col['both'])		
lines(resNoNuggetShape[Dx,],
		resNoNuggetShape[D,],
		type='o', col=col['shape'])		
abline(v=themodel['range']/xres(myraster))


# now with added noise
fracNugget = 1/4
nuggetSd = sqrt(themodel['variance']*fracNugget)
thedf$yNoise = rnorm(nrow(thedf),
		thedf$y,nuggetSd)


source("../R/loglikGmrf.R")

Sar2 = exp(seq(log(0.075),log(0.2),len=12))
Snugget = seq(0.1, 0.3, by=0.02)


if(FALSE) {

resTest = loglikGmrf(
				oneminusar=Sar2Long,
				propNugget = Snugget,
				Yvec=thedf$yNoise,
				Xmat=Xmat,
				shape=maternShape,
				NN=theNN,mc.cores=4,
				adjustEdges=TRUE
		)

		par(mfrow=c(2,2))
for(D in c('logL.ml','logDetVar','rpr','sigmasq.ml')) {		
plot(resTest["rangeInCells",1,],resTest[D,1,],
		ylim=quantile(resTest[D,,],prob=c(0.25,1)),
		ylab=D,type='o', 
lwd=4)		
for(D2 in 2:dim(resTest)[2]) {
	lines(resTest["rangeInCells",D2,],resTest[D,D2,],
		col=D2,type='o')		
}
abline(v=themodel['range']/xres(myraster))
legend(
		'bottomright',legend=Snugget,fill=1:length(Snugget))
}
}

if(FALSE){
oneminusar=0.12
propNugget = 0.2
Yvec=thedf$yNoise
Xmat=Xmat
shape=maternShape
NN=theNN
mc.cores=4
adjustEdges=TRUE
adjustParam=TRUE
adjustShape=FALSE
adjustMarginalVariance=FALSE
}
		
resVanilla = loglikGmrf(
		oneminusar=Sar2 ,
		propNugget = Snugget ,
		Yvec=thedf$yNoise,
		Xmat=Xmat,
		shape=maternShape,
		NN=theNN,mc.cores=4
)


resEdge = loglikGmrf(
		oneminusar=Sar2,		
		propNugget = Snugget,
		Yvec=thedf$yNoise,
		Xmat=Xmat,
		shape=maternShape,
		NN=theNN,mc.cores=4,
		adjustEdges=TRUE,	
		adjustParam=FALSE
)


resAdj = loglikGmrf(
		oneminusar=Sar2,		propNugget = Snugget,
		Yvec=thedf$yNoise,
		Xmat=Xmat,
		shape=maternShape,
		NN=theNN,mc.cores=4,
		adjustEdges=TRUE,	
		adjustParam=TRUE,
		adjustMarginalVariance=TRUE
)


par(mfrow=c(2,2))
D='logL.ml'
Dx='oneminusar'

scale = sqrt(8*themodel['shape'])*xres(myraster)/
		themodel["range"]
a = (scale^2 + 4) 
oneminusar = 1-4/a
rangeInCells = themodel['range']/xres(myraster)
forabline = c(oneminusar, rangeInCells,
	 1/(1+1/fracNugget)		)
names(forabline) = c('oneminusar', 'rangeInCells','propNugget')


colType  = c('resVanilla'='black','resEdge'='blue','resAdj'='red')
for(Dx in c('oneminusar','rangeInCells')) {
	plot(range(resEdge[Dx,,]),c(-15,2) + max(resEdge[D,,]),
			xlab=Dx,ylab=D,type='n')
	for(Dtype in names(colType)) {
		
		res = get(Dtype)
		maxOverNugget = apply(res[D,,],2,which.max)
		lines(diag(res[Dx,maxOverNugget,]),
				diag(res[D,maxOverNugget,]),
col=colType[Dtype],type='o')

}
abline(v=forabline[Dx])
}


Dx='propNugget'

plot(range(resEdge[Dx,,]),quantile(resEdge[D,,],prob=c(0.6,1)),
		xlab=Dx,ylab=D,type='n')
for(Dtype in names(colType)) {
	
	res = get(Dtype)
	maxOverAr = apply(resEdge[D,,],1,which.max)
	lines(diag(res[Dx,,maxOverAr]),
			diag(res[D, ,maxOverAr]),
			col=colType[Dtype],type='o')
	
}
abline(v=forabline[Dx])

D = 'sigmasq.ml'
plot(range(resEdge[Dx,,]),quantile(resEdge[D,,],prob=c(0.2,0.9)),
		xlab=Dx,ylab=D,type='n')
for(Dtype in names(colType)) {
	
	res = get(Dtype)
	maxOverAr = apply(resEdge[D,,],1,which.max)
	lines(diag(res[Dx,,maxOverAr]),
			diag(res[D, ,maxOverAr]),
			col=colType[Dtype],type='o')
	
}

abline(v=forabline[Dx])
abline(h=themodel['variance'])
legend('topright',fill=colType,legend=names(colType))



thesummary = summaryGmrfFit(resAdj) 

themle = thesummary$ml[,'mle']



myraster = raster(extent(0,8000,0,6000), ncols=40,nrows=30)
myrasterBig = extend(myraster, 
		extend(extent(myraster),6*xres(myraster))
)

theNN = NNmat(myraster)
themodel = c(range=5*xres(myraster),shape=2,variance=900)

theU = RFsimulate(myrasterBig,model=themodel)
theU = raster::crop(theU, myraster)

thecov = myraster
values(thecov) = c(rep(1,ncell(thecov)/2),
		rep(4,ncell(thecov)/2))
beta.x=5
theY = theU + beta.x*thecov

theData = brick(theY,thecov)
names(theData) = c("y","x")

fracNugget = 1/4
nuggetSd = sqrt(themodel['variance']*fracNugget)

theData = brick(theY,thecov)
names(theData) = c("y","x")
yNoise = raster(theData)
values(yNoise) = rnorm(ncell(yNoise),
		values(theData[['y']]),nuggetSd)
names(yNoise) = 'yNoise'
theData2 = stack(theData, yNoise)
temp = lgmrfm(theData2, formula = yNoise ~ x,
		oneminusar=exp(seq(log(0.01), log(0.1),len=12)), 
		nugget=seq(0.05, 0.4, len=16),
		NN=theNN,adjustEdges=TRUE,mc.cores=4)
tempV = lgmrfm(theData2, formula = yNoise ~ x,
		oneminusar=exp(seq(log(0.01), log(0.1),len=12)), 
		nugget=seq(0.05, 0.4, len=16),
		NN=theNN,adjustEdges=TRUE,mc.cores=4,
		adjustMarginalVariance=TRUE)

tempA = lgmrfm(theData2, formula = yNoise ~ x,
		oneminusar=exp(seq(log(0.01), log(0.1),len=24)), 
		nugget=seq(0.05, 0.4, len=16),
		NN=theNN,adjustEdges=TRUE,adjustShape=TRUE,
		adjustParam=TRUE,adjustMarginalVariance=TRUE,
		mc.cores=4)
plotLgmrf(temp)
plotLgmrf(tempA)
plotLgmrf(tempV)

