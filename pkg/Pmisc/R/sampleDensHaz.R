#' Sample hazard functions
#' 
#' @description Take posterior samples of hazard funnctions from an INLA fit
#'
#' @param fit result from inla
#' @param x sequence of times
#' @param n number of samples
#' @param covariates values of covariates to use for computing the hazard
#' @param scale scale the hazard (12 converts monthly hazard to yearly)

#' @return matrix or array
#' 
#' @export
sampleDensHaz = function(fit, x,n=1, covariates = c('(Intercept)'=1), scale=1) {
	if(is.vector(covariates)) covariates = t(covariates)
	if(is.null(rownames(covariates)) | 
			all(rownames(covariates) == as.character(1:nrow(covariates)))
		) {
		Srownames = NULL
		for(D in colnames(covariates)[colnames(covariates) != '(Intercept)' ] ) {
			Srownames = paste(Srownames, D, covariates[,D], sep='')
		}	
		getRid = duplicated(Srownames)
		if(any(getRid))
			covariates = covariates[!getRid,,drop=FALSE]
		rownames(covariates) = Srownames[!getRid]
	}
	
	xScale = x/scale
 
	result = replicate(n, sampleDensHazSingle(fit=fit, x=xScale, covariates=covariates))
	if(n < 10)
		dimnames(result)[[length(dim(result))]] = paste("sample", 1:n, sep='')

result[,,'dens',] = result[,,'dens',]/scale

result = drop(result)

result
}

sampleDensHazSingle = function(fit, x, covariates ) {
	samp = INLA::inla.posterior.sample(1,fit)[[1]]
	shape = samp$hyperpar[grep('alpha parameter for ', 
     names(samp$hyperpar))]
	betaIndex = which(fit$misc$configs$contents$length==1)
	betaPars = fit$misc$configs$contents$start[betaIndex]
	beta = samp$latent[betaPars]
	names(beta) = fit$misc$configs$contents$tag[betaIndex]

	if(is.data.frame(covariates)) covariates = as.matrix(covariates)
	if(!is.matrix(covariates)){
		covariates = t(as.matrix(covariates))
	}
	
	covNames = intersect(colnames(covariates), names(beta))
	lambda = drop(exp(-tcrossprod(beta[covNames], 
       covariates[,covNames,drop=FALSE])))
	
	res = abind::abind(
			dens = outer(x, lambda, function(a,b){stats::dweibull(a,shape=shape, scale=b)}),
			cumhaz = outer(x, lambda, '/')^shape,
			along=3
	)
  if(length(x) < 10) dimnames(res)[[1]] = as.character(x)
	res
}



marginalDensHazSingle = function(shape,x){

	cbind(
			haz = sort(as.numeric(
							exp(log(shape) + (shape-1)*log(x)))
			), 
			cumHaz = sort(as.numeric(exp(shape*log(x))))
	)
}

marginalDensHaz = function(fit, x){
	
	if(length(x)==1) x = c(0, x)
	if(length(x)==2){
		x = seq(min(x), max(x), len=200)
	}
	shape = fit$summary.hyperpar[
			grep('alpha parameter for weibull', rownames(fit$summary.hyperpar)), 
			grep("quant$", colnames(fit$summary.hyperpar))]
	
	densHaz = simplify2array(
			mapply(marginalDensHazSingle, x=x, MoreArgs = list(shape=unlist(shape)), SIMPLIFY=FALSE)
	)
	dimnames(densHaz)[[1]] = colnames(shape)
	densHaz		
			
}
