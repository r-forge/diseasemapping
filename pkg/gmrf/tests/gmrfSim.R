library(geostatsp)

myraster = raster(extent(0,8000,0,6000), ncols=40,nrows=30)
theNN = NNmat(myraster)
themodel = c(range=4*xres(myraster),shape=2,var=900)
if(T){
	theU = RFsimulate(myraster,model=themodel)
} else {
	theVar = maternGmrfPrec(theNN, 
		param=c(themodel,cellSize=xres(myraster)),
		adjustEdges=TRUE)



	theChol = chol(theVar)
	theU = as.vector(solve(theChol,rnorm(ncell(myraster))))
}


if(FALSE) {
	
	paramH=c(range=4*xres(myraster),
			shape=2,var=1,cellSize=xres(myraster))
	varMat = matern(myraster,param=paramH)
	
	precEdge = maternGmrfPrec(theNN,param=paramH,
			adjustEdges=TRUE,adjustParam=FALSE)

	precAdj = maternGmrfPrec(theNN,param=paramH,
			adjustEdges=FALSE,adjustParam=TRUE,adjustShape=FALSE)
	
	precBoth = maternGmrfPrec(theNN,param=paramH,
			adjustEdges=TRUE,adjustParam=TRUE,adjustShape=FALSE)
	
	
	prodEdge = precEdge %*% varMat
	prodAdj = precAdj %*% varMat
	prodBoth = precBoth %*% varMat

	
	par(mfrow=c(3,2))
	
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
theY = theU + thecov

theData = brick(theY,thecov)
names(theData) = c("y","x")

thedf = as.data.frame(theData)

Xmat=cbind(inter=1,x=thedf$x)
maternShape=as.numeric(themodel['shape'])

Sar = exp(seq(log(0.075),log(0.25),len=24))

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

		
		
		col = c(vanilla = 'black', adj='green',edge='blue',both='red')
		
		
		scale = sqrt(8*themodel['shape'])*xres(myraster)/
				themodel["range"]
		a = (scale^2 + 4) 
		oneminusar = 1-4/a
		
		par(mfrow=c(2,2))		
		for(D in c('logL.ml','logL.reml')) { 
		
		plot(resNoNuggetEdge['oneminusar',],
				resNoNuggetEdge[D,],
				xlab='1-ar',ylab=D,
				ylim=c(median(resNoNuggetVanilla[D,]),
						max(resNoNuggetBoth[D,])),
				xlim = oneminusar*c(0.25, 1.25),
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
		

		abline(v= oneminusar)

plot(resNoNuggetEdge['rangeInCells',]*xres(myraster),
		resNoNuggetEdge[D,],
		ylim=c(median(resNoNuggetVanilla[D,]),
				max(resNoNuggetBoth[D,])),
		xlim = themodel['range']*c(0.75, 2),
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
abline(v=themodel['range'])
 


}
 	
# now with added noise
nuggetSd = 10
thedf$yNoise = rnorm(nrow(thedf),
		thedf$y,nuggetSd)



source("../R/loglikGmrf.R")

Sar2 = exp(seq(log(0.025),log(0.25),len=4))
Snugget = seq(0,0.9,len=5)
Sar2Long= exp(seq(log(min(Sar2)), log(max(Sar2)),
				len=length(Sar2)*2))
SnuggetLong = seq( min(Snugget) ,  max(Snugget) ,
				len=length(Snugget)*2)



resVanilla = loglikGmrf(
		oneminusar=Sar2Long,
		propNugget = SnuggetLong,
		Yvec=thedf$yNoise,
		Xmat=Xmat,
		shape=maternShape,
		NN=theNN,mc.cores=4
)

par(mfrow=c(2,2))

plot(resVanilla['oneminusar',1,],
		resVanilla['logL',1,])


plot(resVanilla['oneminusar',1,],
		apply(resVanilla['logL',,],2,max))

plot(resVanilla['propNugget',,1],
		apply(resVanilla['logL',,],1,max))

image(resVanilla['propNugget',,1], resVanilla['oneminusar',1,], 
		resVanilla['logL',,],col=terrain.colors(12))

resAdj = loglikGmrf(
		oneminusar=Sar2Long,		
		propNugget = SnuggetLong,
		Yvec=thedf$yNoise,
		Xmat=Xmat,
		shape=maternShape,
		NN=theNN,mc.cores=4,
		adjustParam=TRUE,adjustEdges=FALSE
)

resEdge = loglikGmrf(
		oneminusar=Sar2,		
		propNugget = Snugget,
		Yvec=thedf$yNoise,
		Xmat=Xmat,
		shape=maternShape,
		NN=theNN,mc.cores=4,
		adjustEdges=TRUE,	adjustParam=FALSE
)


resBoth = loglikGmrf(
		oneminusar=Sar2,		propNugget = Snugget,
		Yvec=thedf$yNoise,
		Xmat=Xmat,
		shape=maternShape,
		NN=theNN,mc.cores=4,
		adjustEdges=TRUE,	adjustParam=TRUE
)




thesummary = summaryGmrfFit(res1,theArgs)$summary
thesummary

plot(res1['rangeInCells',],res1['logL',],type='o',log='x',
		ylim=c(-20,0)+max(res1['logL',]))
Nar = 6
text(res1['rangeInCells',
				round(seq(1,ncol(res1),len=Nar))],
		min(res1['logL',]),
		signif(res1['ar',round(seq(1,ncol(res1),len=Nar))],3)
)
abline(v=themodel['range']/xres(myraster),col='blue')
abline(v=thesummary['rangeInCells',
				c('mle','q0.025','q0.975')],
		col=c('red','orange','orange'))

