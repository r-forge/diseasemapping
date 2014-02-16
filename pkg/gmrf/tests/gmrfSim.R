library(geostatsp)

myraster = raster(extent(0,8000,0,6000), ncols=40,nrows=30)
theNN = NNmat(myraster)
themodel = c(range=12*xres(myraster),shape=2,var=900)
if(TRUE){
	theU = RFsimulate(myraster,model=themodel)
} else {
	theVar = maternGmrfPrec(theNN, 
		param=c(themodel,cellSize=xres(myraster)),
		adjustEdges=FALSE)

	theChol = chol(theVar)
	theU = as.vector(solve(theChol,rnorm(ncell(myraster))))
}
thecov = myraster
values(thecov) = c(rep(1,ncell(thecov)/2),
		rep(4,ncell(thecov)/2))
theY = theU + thecov

theY = brick(theY,thecov)

thedf = as.data.frame(theY)

theArgs = list(
		Yvec=thedf[,1],
		Xmat=cbind(inter=1,x=thedf[,2]),
		NN=theNN,
		maternShape=as.numeric(themodel['shape'])
)

SrangeInCells=themodel['range']*
		exp(seq(log(0.5),log(2),len=12))/
		xres(myraster) 

res1 = parallel::mcmapply(loglikGmrf,
		rangeInCells=SrangeInCells, 
		MoreArgs=theArgs,mc.cores=4)


thesummary = summaryGmrfFit(res1,theArgs)$summary

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
thesummary


# with edge correction

theArgs$adjustEdges=TRUE
res2 = parallel::mcmapply(loglikGmrf,
		rangeInCells=SrangeInCells, 
		MoreArgs=theArgs,mc.cores=4)
lines(res2['rangeInCells',],res2['logL',],type='o',
		col='green')
abline(v=thesummary['rangeInCells',
				c('mle','q0.025','q0.975')],
		col=c('green','yellow','yellow'))


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

