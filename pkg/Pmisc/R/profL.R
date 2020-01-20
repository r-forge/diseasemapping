#' Likelihood surfaces
#' 
#' @description Evalues multivariable profile likelihoods
#'
#' @param fit model fit from, for example, glm
#' @param se number of standard errors from MLE of each parameter
#' @return List which includes an array of likelihoods
#' @examples
#'     # from glm  help file
#'     counts <- c(18,17,15,20,10,20,25,13,12)
#'     outcome <- gl(3,1,9)
#'     treatment <- gl(3,3)
#'     glm.D93 <- glm(counts ~ outcome, family = poisson())
#' x = likSurface(glm.D93)
#' image(x$parameters[,1], x$parameters[,2], apply(x$logLik, c(2,3), min))
#' 
#' @export
likSurface = function(fit, se=seq(0,3,len=10)){
	
	parameters = outerParameters(fit, se)
	
	family = as.character(fit$family)[1]

	if(family=='poisson') {
		likFun = likFunPoisson
	} else if(family=='binomial'){
		likFun = likFunBinom
	} else {
		stop("family", family, 'not available')
	}
		

	yCol = all.vars(fit$terms)[1]


	likVec = apply(parameters$long, 1, likFun, 
			y=fit$model[,yCol], 
			covariates=stats::model.matrix(fit$formula, fit$model))
	
	likArray = array(likVec, rep(nrow(parameters$short), length(fit$coef)))
	
	likBreaks = likelihoodBreaks(
			max(likArray),  df=length(fit$coef))
	
	result = c(list(
					logLik = likArray,
					parameters = parameters$short,
					mle = fit$coef),
			infMat(fit, y=fit$model[[yCol]], 
				covariates=stats::model.matrix(
					fit$formula, fit$model), 
				likFun), 
			likBreaks)

	if(ncol(parameters$short)>1) {
	for(D in seq(1, ncol(parameters$short)-1)) {
		result$profileBreaks[[paste("df", D, sep='')]] =
				max(result$logLik) - stats::qchisq(result$breaks, df=D)/2
	}
	}
	result	
}


likFunBinom = function(param, y, covariates){
	
	x = param
	xNotIntercept = x[grep("^\\(Intercept\\)$", names(x), invert=TRUE)]
	logitmu = x['(Intercept)'] + 
			as.matrix(covariates[,names(xNotIntercept)]) %*% xNotIntercept
	
	mu = exp(logitmu)/(1+exp(logitmu))
	
	sum(stats::dbinom(y[,1], y[,1]+y[,2], mu, log=TRUE))  
}

likFunPoisson = function(param, y, covariates){
	
	x = param
	xNotIntercept = x[grep("^\\(Intercept\\)$", names(x), invert=TRUE)]
	logmu = x['(Intercept)'] + 
			as.matrix(covariates[,names(xNotIntercept)]) %*% xNotIntercept
	
	offsetCols = grep("^offset\\(", colnames(covariates), value=TRUE)
	for(D in offsetCols){
		logmu = logmu+ covariates[,offsetCols]
	}
	
	mu = exp(logmu)
	
	sum(stats::dpois(y, mu, log=TRUE))  
}


infMat = function(fit, y, covariates, likFun) {
		hessian =  numDeriv::hessian(
				likFun,
				fit$coef,
				y=y, covariates=covariates
				)
				
		dimnames(hessian) = list(
				names(fit$coefficients),
				names(fit$coefficients)
				)		
				
		result = list(
			hessian = hessian, 
			infMat = -hessian
		)

		result$varParam = solve(result$infMat)
		result$se = sqrt(diag(result$varParam))
				
		result
}

likelihoodBreaks = function(
		logLik=0,
		prob = c(0, 0.5, 0.8, 0.95, 0.99, 0.999),
		df=1
) {

	prob = sort(prob, decreasing=TRUE)
	Squant = max(logLik) - stats::qchisq(prob, df=df)/2
	Scol = RColorBrewer::brewer.pal(length(Squant)-1, 'YlGnBu')
	
	list(quant=Squant, col=Scol, breaks=prob, df=df)
	
}

outerParameters = function(
		x,	se
) {
	
	Sse = sort(unique(c(-se, se)))
	parMat = outer(summary(x)$coef[,'Std. Error'], Sse)
	parMat = t(parMat + x$coef)
	parMatX = do.call(expand.grid, as.data.frame(parMat))
	
	list(short=parMat, long=as.matrix(parMatX))
}

