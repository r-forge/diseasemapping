library(geostatsp)

myraster = squareRaster(raster(extent(0,8000,0,6000), ncols=80,nrows=60))




myQ = maternGmrfPrec(myraster, 
		param=c(shape=2, oneminusar=0.1,
				conditionalVariance=100),
		adjustEdges='optimalShape')
library(RColorBrewer)
thecol = c('#000000', 
		brewer.pal(ncol(attributes(myQ)$par$emp)-2,'Set2')
	)
matplot(attributes(myQ)$par$emp$x,
		attributes(myQ)$par$emp[,-1],
		type='l',
#		log='y',#ylim=c(500, 1000),xlim=c(0,1000),
		lty=1, lwd=(7:2)/3, col=thecol
				)

legend('topright',lty=1, col=thecol, 
		legend=colnames(attributes(myQ)$par$emp)[-1])


themodel = attributes(myQ)$param$optimalShape
maternShape = attributes(myQ)$param$theo['shape']

theU = brick(myraster, nl=200)
for(D in 1:nlayers(theU))
	theU[[D]] = RFsimulate(myraster,model=themodel)

if(FALSE){
	maternmat = matern(myraster, #param=themodel)
			param=)
#	diag(maternmat)= sum(attributes(myQ)$param$optimalWithNugget[
#					c('variance','nugget')])
	maternchol = chol(maternmat, LDL=FALSE,pivot=FALSE)
	for(D in 1:nlayers(theU))
	values(theU[[D]]) = 
			as.numeric(crossprod(maternchol,
			rnorm(ncol(maternmat))))
}

if(FALSE){
	temp = (myQ %*% matern(
						myraster, #param=themodel)
						param=
								attributes(myQ)$param$optimal))
	par(mfrow=c(1,2))
	hist(diag(temp),breaks=100)
	hist(temp[lower.tri(temp,diag=F)],breaks=100)
}


thecov = myraster
values(thecov) = c(rep(0,ncell(thecov)/2),
		rep(4,ncell(thecov)/2))
names(thecov)='x'
beta.x=5
theY = theU + beta.x*thecov

 
theNN = NNmat(theY)



Sar = exp(seq(log(0.05),log(0.15),len=56))

source("../R/loglikGmrf.R")
source("../R/lgmrfm.R")
source("../R/conditionalGmrf.R")
source("../R/loglikGmrf2.R")


 

 

resNoNugget = 
		loglikGmrf(
		oneminusar=Sar ,
		propNugget = 0 ,
		Yvec=as.data.frame(theY),
		Xmat=cbind(intercept=1,as.data.frame(thecov)),
		shape=maternShape,
		NN=theNN,mc.cores=4
)

plotres = function(res) {

	par(mfrow=c(1,2))
	logL = res['logL.reml',,]
	themax = apply(logL, 1, max)
	logL = logL - array(themax, dim(logL))
	logL = t(logL)

	
	thex = res['range',1,]
	matplot(thex, logL, xlab='range', ylab='logL',
			type='l', lty=1, col='#00000030',
			ylim=c(-4, 0)#,
#			xlim = c(0.8, 1.25)*attributes(myQ)$param$theo['range']
	)
	abline(h=-2,col='yellow')
	
	themax = apply(logL, 2, which.max)
points(
		thex[themax], 
		rep(-2, ncol(logL)), col='#FF000030',
		pch=16) 
		abline(v=attributes(myQ)$param$theo['range'],
				col='red')

thepar=		'xisq.reml'
toplot = 		res[thepar, , themax]
thetrue =attributes(myQ)$param$original['conditionalVariance'] 
hist(toplot, main='', prob=TRUE, 
		xlim=range(c(toplot, thetrue)),
		xlab=thepar)	
abline(v=thetrue,col='red')
}

plotres(resNoNugget)



resNoNuggetEdge = loglikGmrf(
		oneminusar=Sar ,
		propNugget = 0 ,
		Yvec=as.data.frame(theY),
		Xmat=cbind(intercept=1,as.data.frame(thecov)),
		shape=maternShape,
		NN=theNN,mc.cores=4,
		adjustEdges='optimal'
)
plotres(resNoNuggetEdge)

 

if(FALSE){
oneminusar=0.12
propNugget = 0.2
Yvec=thedf$yNoise
Xmat=Xmat
shape=maternShape
NN=theNN
mc.cores=4
adjustEdges=FALSE
}

# now with added noise
fracNugget = 1/2
nuggetSd = sqrt(themodel['variance']*fracNugget)
yNoise = theY
values(yNoise) = values(yNoise) + 
		rnorm(nlayers(yNoise)*ncell(yNoise),
				mean=0, sd=nuggetSd)


nuggetSd^2/
		attributes(myQ)$par$theo['conditionalVariance']

Sar2 =  seq((0.04),(0.2),len=40)
Snugget =  seq((0.025), (0.4),len=100)
 


resVanilla = loglikGmrf(
		oneminusar=Sar2 ,
		propNugget = Snugget ,
		Yvec=as.data.frame(yNoise),
		Xmat=cbind(intercept=1,as.data.frame(thecov)),
		shape=maternShape,
		NN=theNN,mc.cores=4
)
 
bob = function(res) { 
dseq = rev(c(0,0.5, 1,2,4,8,20))

res = res[,,seq(dim(res)[3],1)]

thecol = mapmisc::colourScale(1,
		breaks=c(min(res['logL.ml',,]), 
				max(res['logL.ml',,]) - dseq), 
		col='RdYlGn',style='fixed',rev=TRUE)

plot(range(res['propNugget',,1]), 
		range(res['range',1,]),type='n',
		xlab='tausq/xisq',ylab='range',log='xy')
.filled.contour(res['propNugget',,1], 
		res['range',1,],res['logL.ml',,],
		col=thecol$col,levels=thecol$breaks
		)
		
points(nuggetSd^2/
				attributes(myQ)$par$theo['conditionalVariance'],
		themodel['range'], cex=2)		
return(list(dseq=dseq, col=thecol$col))		
}

bob2 = function(res) {
	
	par(mfrow=c(2,1))
	logL = apply(res['logL.reml',,,], 
			c(1,3), max)
	themax = apply(logL, 1, max)
	logL = logL - array(themax, dim(logL))
	logL = t(logL)
	
	thex = res['range',1,1,]
	
	matplot(thex,
			logL,
			lty=1, type='l',
			col='#00000030',
			xlab='range',ylab='logL',
			ylim=c(-4,0))
	
	
	themax = apply(logL, 2, which.max)
	denstrip::denstrip(
			thex[themax], 
			at=-4, 
			colmax='red',
			ticks = quantile(thex[themax],
					prob=c(0.2, 0.5, 0.8)),
			tcol='orange') 
	abline(v=attributes(myQ)$param$theo['range'],
			col='red')
	abline(h=-2,col='yellow')
	
	logL = apply(res['logL.reml',,,], 
			c(1,2), max)
	themax = apply(logL, 1, max)
	logL = logL - array(themax, dim(logL))
	logL = t(logL)
	
	thex = res['propNugget',1, ,1]
	
	matplot(thex,
			logL,
			lty=1, type='l',
			col='#00000030',
			xlab='tausq/xisq',ylab='logL',
			ylim=c(-4,0))
	
	
	themax = apply(logL, 2, which.max)
	denstrip::denstrip(
			thex[themax], 
			at=-4, 
			colmax='red',
			ticks = quantile(thex[themax],
					prob=c(0.2, 0.5, 0.8)),
			tcol='orange') 
	
	abline(v=nuggetSd^2/attributes(myQ)$param$theo['conditionalVariance'],
			col='red')
	abline(h=-2,col='yellow')

	
	
}

par(mfrow=c(4,4),mar=c(2.2,2.2,0,0))
for(D in 1:prod(par('mfrow')))
	forlegend = bob(resVanilla[,D,,])
mapmisc::legendBreaks("bottomright",
		col=forlegend$col,breaks=forlegend$dseq)

bob2(resVanilla)

resEdge = loglikGmrf(
		oneminusar=Sar2,		
		propNugget = Snugget,
		Yvec=as.data.frame(yNoise),
		Xmat=cbind(intercept=1,as.data.frame(thecov)),
		shape=maternShape,
		NN=theNN,mc.cores=4,
		adjustEdges='optimal'
)
#pdf("/tmp/stuff.pdf")
par(mfrow=c(4,4),mar=c(2.2,2.2,0,0))
for(D in 1:prod(par('mfrow')))
	forlegend = bob(resEdge[,D,,])
mapmisc::legendBreaks("topleft",
		col=forlegend$col,breaks=forlegend$dseq)
#dev.off()
bob2(resEdge)


if(FALSE ) {

resNoNugget  = lgmrfm(
		data=theData,
		formula = y ~ x,
		oneminusar=Sar,
		shape=maternShape,
		mc.cores=4,
		NN = theNN
)
plot(resNoNugget["oneminusar",],resNoNugget["logL.ml",])


resNoNuggetEdge = lgmrfm(
		data=theData,
		formula = y ~ x,
		oneminusar=Sar,
		shape=maternShape,
		mc.cores=4,
		NN = theNN,
		adjustEdges=TRUE
)

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
}