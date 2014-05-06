library(geostatsp)

myraster = squareRaster(raster(extent(0,8000,0,6000), ncols=60,nrows=40))
myrasterBig = extend(myraster, 
		extend(extent(myraster),6*xres(myraster))
)

themodel = c(range=3*xres(myraster),shape=2,variance=900)

theU = RFsimulate(myrasterBig,model=themodel)
theU = raster::crop(theU, myraster)

	
thecov = myraster
values(thecov) = c(rep(1,ncell(thecov)/2),
		rep(4,ncell(thecov)/2))
beta.x=5
theY = theU + beta.x*thecov

theData = brick(theY,thecov)
names(theData) = c("y","x")
theNN = NNmat(theData)


Sar = exp(seq(log(0.04),log(0.16),len=12))


source("../R/loglikGmrf.R")
source("../R/lgmrfm.R")
source("../R/conditionalGmrf.R")


 

resNoNugget  = lgmrfm(
		data=theData,
		formula = y ~ x,
		oneminusar=Sar,
		shape=themodel['shape'],
		mc.cores=2,
		NN = theNN
		)

		

resNoNuggetEdge = lgmrfm(
		data=theData,
		formula = y ~ x,
		oneminusar=Sar,
		shape=themodel['shape'],
		mc.cores=1,
		NN = theNN,
				adjustEdges=TRUE
		)

if(FALSE) {
		resNoNuggetAdj = lgmrfm(
				data=theData,
				formula = y ~ x,
				oneminusar=Sar,
				shape=themodel['shape'],
				mc.cores=2,
				NN = theNN,
				adjustEdges=TRUE,
				adjustMarginalVariance=TRUE
		)
	} else{
		resNoNuggetAdj = NULL
	}
		
		col = c(vanilla = 'black', adj='green',edge='blue',both='red',shape='orange')
		
		
		scale = sqrt(8*themodel['shape'])*xres(myraster)/
				themodel["range"]
		a = (scale^2 + 4) 
		oneminusar = 1-4/a
		
		par(mfrow=c(2,2))		
		D ='logL.ml' 
		
		matplot(
				resNoNuggetEdge$complete['oneminusar',],
				cbind(
						resNoNuggetEdge$complete[D,],
						resNoNugget$complete[D,]),
				xlab='1-ar',ylab=D,
				col=col[c('edge','vanilla')],type='o',
				log='x',lwd=1,pch=16)
		legend('bottom',fill=col[c('vanilla','edge')],
				legend=c('vanilla','edge'))
		abline(v= oneminusar,col='red')
		theci = c('mle'=1,'q0.025'=3,'q0.975'=3)
		abline(v=resNoNuggetEdge$ml['oneminusar',names(theci)],
				lty=theci,col=col['edge'])
		abline(v=resNoNugget$ml['oneminusar',names(theci)],
				lty=theci,col=col['vanilla'])
		
		
		matplot(
				resNoNuggetEdge$complete['rangeInCells',],
				cbind(
						resNoNuggetEdge$complete[D,],
						resNoNugget$complete[D,]),
				xlab='range',ylab=D,
				col=col[c('edge','vanilla')],type='o',
				log='x',lwd=1,pch=16)

abline(v=themodel['range']/xres(theData	),col='red')

abline(v=resNoNuggetEdge$ml['rangeInCells',names(theci)],
		lty=theci,col=col['edge'])
abline(v=resNoNugget$ml['rangeInCells',names(theci)],
		lty=theci,col=col['vanilla'])



matplot(
		resNoNuggetEdge$complete['rangeInCells',],
		cbind(
				resNoNuggetEdge$complete[D,],
				resNoNugget$complete[D,]),
		xlab='range post fit',ylab=D,
		col=col[c('edge','vanilla')], 
		log='x',lwd=1,pch=16,type='n')
 
matlines(
		resNoNuggetEdge$complete['range.postfit',],
		resNoNuggetEdge$complete[D,],
		col=col['edge'],
		lty=2,pch=1,type='o'
)
matlines(
		resNoNugget$complete['range.postfit',],
		resNoNugget$complete[D,],
		col=col['vanilla'],
		lty=2,pch=1,type='o'
)

abline(v=themodel['range']/xres(theData	),col='red')

abline(v=resNoNuggetEdge$ml['range.postfit',names(theci)],
		lty=theci,col=col['edge'])
abline(v=resNoNugget$ml['range.postfit',names(theci)],
		lty=theci,col=col['vanilla'])


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
adjustEdges=FALSE
adjustParam=FALSE
adjustShape=FALSE
adjustMarginalVariance=FALSE
}

# now with added noise
fracNugget = 1/10
nuggetSd = sqrt(themodel['variance']*fracNugget)
thedf = as.data.frame(theData)
thedf$intercept = 1
Xmat = as.matrix(thedf[,c("intercept","x")])
thedf$yNoise = rnorm(nrow(thedf),
		thedf$y,nuggetSd)


Sar2 = exp(seq(log(0.01),log(0.05),len=12))
Snugget = seq(5, 50,len=25)
maternShape = themodel['shape']


resVanilla = loglikGmrf(
		oneminusar=Sar2 ,
		propNugget = Snugget ,
		Yvec=thedf$yNoise,
		Xmat=Xmat,
		shape=maternShape,
		NN=theNN,mc.cores=4
)
 
image(resVanilla['propNugget',,1], resVanilla['oneminusar',1,],resVanilla['logL.ml',,])


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

fracNugget = 1/10
nuggetSd = sqrt(themodel['variance']*fracNugget)

theData = brick(theY,thecov)
names(theData) = c("y","x")
yNoise = raster(theData)
values(yNoise) = rnorm(ncell(yNoise),
		values(theData[['y']]),nuggetSd)
names(yNoise) = 'yNoise'
theData2 = stack(theData, yNoise)
if(FALSE){
	data=theData2
	formula = yNoise ~x
	oneminusar=0.2
	oneminusar=Sar
	nugget=0.2
	nugget=Snugget
	shape=themodel['shape']
	NN=theNN;adjustEdges=TRUE
	
	dataOrig=data
	data = as.data.frame(data)
	Yvec = data[,
			as.character(attributes(terms(formula))$variables)[2]
	]
	Xmat = model.matrix(formula,data=data)
	
	Yvec=Yvec
	Xmat=Xmat 
	NN=NN
	propNugget=nugget  
	shape=shape
	
	param=mleparam
	template=dataOrig
	
	
}


source("../R/lgmrfm.R")
source("../R/conditionalGmrf.R")
Snugget = exp(seq(log(0.01), log(0.2), len=16))
Sar = exp(seq(log(0.05), log(0.25),len=25))
temp = lgmrfm(theData2, formula = yNoise ~ x,
		oneminusar=Sar, 
		nugget=Snugget,shape=themodel['shape'],
		NN=theNN,adjustEdges=FALSE,mc.cores=4)
tempV = lgmrfm(theData2, formula = yNoise ~ x,
		oneminusar=Sar, 
		nugget=Snugget,shape=themodel['shape'],
		NN=theNN,adjustEdges=TRUE,mc.cores=4)

tempA = lgmrfm(theData2, formula = yNoise ~ x,
		oneminusar=Sar, 
		nugget=Snugget,shape=themodel['shape'],
		NN=theNN,adjustEdges=TRUE,adjustShape=FALSE,
		adjustParam=TRUE,adjustMarginalVariance=TRUE,
		mc.cores=4)

plotLgmrf(temp)

par(mfrow=c(2,2))
plot(theU)
plot(temp$predict[['random']])
plot(temp$predict[['krigeSd']])
plot(temp$predict[['fixed']])
