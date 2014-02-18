library(geostatsp)

myraster = raster(extent(0,8000,0,6000), ncols=40,nrows=30)
theNN = NNmat(myraster)
themodel = c(range=6*xres(myraster),shape=2,var=900)
if(T){
	theU = RFsimulate(myraster,model=themodel)
} else {
	theVar = maternGmrfPrec(theNN, 
		param=c(themodel,cellSize=xres(myraster)),
		adjustEdges=TRUE)



	theChol = chol(theVar)
	theU = as.vector(solve(theChol,rnorm(ncell(myraster))))
}
if(FALSE){
	theV = matern(myraster, param=themodel)
	theC = chol(theV)
	theU = theC %*% rnorm(ncell(myraster))
}

thecov = myraster
values(thecov) = c(rep(1,ncell(thecov)/2),
		rep(4,ncell(thecov)/2))
theY = theU + thecov

theY = brick(theY,thecov)
names(theY) = c("y","x")

thedf = as.data.frame(theY)

theArgs = list(
		Yvec=thedf$y,
		Xmat=cbind(inter=1,x=thedf$x),
		NN=theNN,
		maternShape=as.numeric(themodel['shape'])
)

SrangeInCells=themodel['range']*
		seq(0.75,1.5,len=24)/
		xres(myraster) 

source("../R/loglikGmrf.R")

res1 = parallel::mcmapply(loglikGmrf,
		rangeInCells=SrangeInCells, 
		MoreArgs=theArgs,mc.cores=4)


thesummary = summaryGmrfFit(res1,theArgs)$summary

# with edge correction
theArgs2 = theArgs
theArgs2$adjustEdges=TRUE

res2 = parallel::mcmapply(loglikGmrf,
		rangeInCells=SrangeInCells, 
		MoreArgs=theArgs2,mc.cores=4)

thesummary2 = summaryGmrfFit(res2,theArgs)$summary


plot(res1['rangeInCells',],res1['logL',],type='o',log='x',
		ylim=c(-6,0)+max(res1['logL',]),
		xlab='range (in cells)',ylab='logL')
Nar = 6
text(res1['rangeInCells',
				round(seq(1,ncol(res1),len=Nar))],
		min(res1['logL',]),
		signif(res1['ar',round(seq(1,ncol(res1),len=Nar))],3)
)
abline(v=themodel['range']/xres(myraster),col='blue')
abline(v=thesummary['rangeInCells',
				c('mle','q0.025','q0.975')],
		col=c('red','orange','orange'),lwd=2)



lines(res2['rangeInCells',],
		res2['logL',]-max(res2['logL',])+
				thesummary['logL','mle'],
		type='o',
		col='green')
abline(v=thesummary2['rangeInCells',
				c('mle','q0.025','q0.975')],
		col=c('green','yellow','yellow'))

thesummary[c('rangeInCells','sigmasq'),]
thesummary2[c('rangeInCells','sigmasq'),]



# now with added noise

theArgsN=theArgs
theArgsN$Yvec = rnorm(length(theArgsN$Yvec),
			theArgs$Yvec, 100)
theArgsN$adjustEdges=FALSE
theArgsN$propNugget=seq(0,0.5,len=10)

res3 = parallel::mcmapply(loglikGmrf,
		rangeInCells=exp(seq(log(4),log(20),len=24)), 
		MoreArgs=theArgs,mc.cores=4)


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

