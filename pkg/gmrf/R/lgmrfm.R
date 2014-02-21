lgmrfm = function(data,formula,covariates=NULL,
		shape=1,nugget=seq(0,1,len=5),
		oneminusar=seq(0.01, 0.4,len=4), 
		rangeInCells=NULL,
		NN=NNmat(data),mc.cores=1,
		...) {

	if(!is.null(covariates)){
		covariates = stackRasterList(covariates,template=data)
		data = stack(data,covariates)
	}
	data = as.data.frame(data)
	Yvec = data[,
			as.character(attributes(terms(formula))$variables)[2]
	]
	Xmat = model.matrix(formula,data=data)
		
	thel = loglikGmrf(oneminusar=oneminusar, rangeInCells=NULL,
			Yvec=Yvec,Xmat=Xmat,
			NN=NN,propNugget=nugget,
			shape=shape,mc.cores=mc.cores,...)
	
	thesummary = summaryGmrfFit(thel)
	thesummary$complete = thel
	
	return(thesummary)
}

plotLgmrf = function(x,reml=FALSE){
	oldpar = par()
	par(mfrow=c(1,2))
	plotcols = c('mle','q0.025','q0.975')
	
	 plot(x$profL$rangeInCells[[c('ml','reml')[1+reml] ]],type='o')
	 abline(v=x[[c('ml','reml')[1+reml]]][
					 'rangeInCells',plotcols
					 ],col=c('red','orange','orange'))
	 
	 plot(x$profL$propNugget[,c(1,2+reml)],type='o')
	 abline(v=x[[c('ml','reml')[1+reml]]][
					 'propNugget',plotcols
			 ],col=c('red','orange','orange'))
	 par(oldpar)
}








