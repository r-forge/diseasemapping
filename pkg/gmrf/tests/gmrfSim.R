
if(Sys.info()['nodename'] == 'darjeeling') {
	ncores = 20
} else if(Sys.info()['nodename'] == 'mud'){
	ncores = 4
} else {
	ncores = 1
}


library(geostatsp)

myraster = squareRaster(raster(extent(0,8000,0,6000), ncols=60,nrows=45))



myQ = maternGmrfPrec(myraster, 
		param=c(shape=2, oneminusar=0.1,
				conditionalVariance=100),
		adjustEdges=FALSE)
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


themodel = attributes(myQ)$param$optimal #Shape
maternShape = attributes(myQ)$param$theo['shape']

theU = RFsimulate(myraster, model=themodel, n=250)
 

if(FALSE){
	stuff='optimal'
	temp = (maternGmrfPrec(myraster, 
						param=c(shape=2, oneminusar=0.1,
								conditionalVariance=100),
						adjustEdges=stuff) %*% matern(
						myraster, 
						param= attributes(myQ)$param[[stuff]]))
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


#source("../R/loglikGmrf.R")
#source("../R/lgmrfm.R")
#source("../R/conditionalGmrf.R")
source("../R/loglikGmrf2.R")

 
Sar = exp(seq(log(0.05),log(0.15),len=40))

resNoNugget = 
		loglikGmrf(
		oneminusar=Sar ,
		propNugget = 0 ,
		Yvec=as.data.frame(theY),
		Xmat=cbind(intercept=1,as.data.frame(thecov)),
		shape=maternShape,
		NN=theNN,mc.cores=ncores
)

plotres = function(res) {

	par(mfrow=c(1,2))
	logL = res['logL.reml',,]
	themax = apply(logL, 1, max)
	logL = logL - array(themax, dim(logL))
	logL = t(logL)

	
	thex = res['oneminusar',1,]
	matplot(thex, logL, xlab='oneminusar', ylab='logL',
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
		abline(v=attributes(myQ)$param$theo['oneminusar'],
				col='red')

thepar=		'xisq.reml'
toplot = 		res[thepar, , themax]
thetrue =attributes(myQ)$param$original['conditionalVariance'] 
hist(toplot, main='', prob=TRUE, 
		xlim=range(c(toplot, thetrue)),
		xlab=thepar)	
abline(v=thetrue,col='red')
}

pdf("resNoNugget.pdf")
plotres(resNoNugget)
dev.off()


resNoNuggetEdge = loglikGmrf(
		oneminusar=Sar ,
		propNugget = 0 ,
		Yvec=as.data.frame(theY),
		Xmat=cbind(intercept=1,as.data.frame(thecov)),
		shape=maternShape,
		NN=theNN,mc.cores=ncores,
		adjustEdges='optimal'
)
pdf("resNoNuggetEdge.pdf")
plotres(resNoNuggetEdge)
dev.off()
 

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
 
nuggetSd =  4  
nuggetSd^2/
		attributes(myQ)$par$theo['conditionalVariance']
yNoise = theY
values(yNoise) = values(yNoise) + 
		rnorm(nlayers(yNoise)*ncell(yNoise),
				mean=0, sd=nuggetSd)




Sar2 =  seq((0.05),(0.25),by=0.005)
Snugget =  seq((0.01), (0.3),by=0.005)

resVanilla = loglikGmrf(
		oneminusar=Sar2 ,
		propNugget = Snugget ,
		Yvec=as.data.frame(yNoise),
		Xmat=cbind(intercept=1,as.data.frame(thecov)),
		shape=maternShape,
		NN=theNN,mc.cores=ncores
)
 
bob = function(res) { 
dseq = rev(c(0,0.5, 1,2,4,8,20))

#res = res[,,seq(dim(res)[3],1)]

thecol = mapmisc::colourScale(1,
		breaks=c(min(res['logL.ml',,]), 
				max(res['logL.ml',,]) - dseq), 
		col='RdYlGn',style='fixed',rev=TRUE)

plot(range(res['propNugget',,1]), 
		range(res['oneminusar',1,]),type='n',
		xlab='tausq/xisq',ylab='oneminusar',log='x')
.filled.contour(res['propNugget',,1], 
		res['oneminusar',1,],res['logL.ml',,],
		col=thecol$col,levels=thecol$breaks
		)
		
points(nuggetSd^2/
				attributes(myQ)$par$theo['conditionalVariance'],
		attributes(myQ)$param$theo['oneminusar'], cex=2)		
return(list(dseq=dseq, col=thecol$col))		
}

bob2 = function(res) {
	
	par(mfrow=c(2,1),mar=c(3,2,0.1, 0.1), mgp = c(2,1,0))
	logL = apply(res['logL.reml',,,], 
			c(1,3), max)
	themax = apply(logL, 1, max)
	logL = logL - array(themax, dim(logL))
	logL = t(logL)
	
	thex = res['oneminusar',1,1,]
	
	matplot(thex,
			logL,
			lty=1, type='l',
			col='#00000030',
			xlab='oneminusar',ylab='logL',
			ylim=c(-4,0))
	
	
	themax = apply(logL, 2, which.max)
	denstrip::denstrip(
			thex[themax], 
			at=-4, 
			colmax='red',
			ticks = quantile(thex[themax],
					prob=c(0.2, 0.5, 0.8)),
			tcol='orange') 
	abline(v=attributes(myQ)$param$theo['oneminusar'],
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
			ylim=c(-4,0), log='x')
	
	
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

pdf("resVanilla1.pdf")
par(mfrow=c(4,4),mar=c(2.2,2.2,0,0), oma=c(1,1,0,0))
for(D in 1:prod(par('mfrow')))
	forlegend = bob(resVanilla[,D,,])
mapmisc::legendBreaks("bottomright",
		col=forlegend$col,breaks=forlegend$dseq)
mtext('oneminusar', side=2, outer=T)
mtext('tausq/xisq', side=1, outer=T)
dev.off()

pdf("resVanilla2.pdf")
bob2(resVanilla)
dev.off()

Sar2 =  seq((0.025),(0.175),by=0.005)
Snugget =  seq((0.025), (0.5),by=0.005)
resEdge = loglikGmrf(
		oneminusar=Sar2,		
		propNugget = Snugget,
		Yvec=as.data.frame(yNoise),
		Xmat=cbind(intercept=1,as.data.frame(thecov)),
		shape=maternShape,
		NN=theNN,mc.cores=ncores,
		adjustEdges='optimal'
)
pdf("resEdge1.pdf")
par(mfrow=c(4,4),mar=c(2.2,2.2,0,0), oma=c(1,1,0,0))
for(D in 1:prod(par('mfrow')))
	forlegend = bob(resEdge[,D,,])
mapmisc::legendBreaks("topleft",
		col=forlegend$col,breaks=forlegend$dseq)
mtext('oneminusar', side=2, outer=T)
mtext('tausq/xisq', side=1, outer=T)
dev.off()

pdf("resEdge2.pdf")
bob2(resEdge)
dev.off()


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