

informationLgm = function(fit, ...) {
	nonLinearParams = c('boxcox','shape','nugget','variance',
			'anisoAngleDegrees','anisoRatio','range')
	
	reEstimate = rownames(fit$summary)[
			fit$summary[,"Estimated"]
	]
	reEstimate = gsub("sdNugget", "nugget", reEstimate)
	reEstimate = gsub("sdSpatial", "variance", reEstimate)
	reEstimate = gsub("range \\(km\\)", "range", reEstimate)
	reEstimate = intersect(reEstimate, nonLinearParams)
	
	baseParam = fit$param[reEstimate]
	moreParams = fit$param[
			!names(fit$param) %in% reEstimate &
					names(fit$param) %in% nonLinearParams]
	
	hess = numDeriv::hessian(geostatsp::loglikLgm, baseParam,
			data=fit$data,trend=fit$model$trend,
			reml=fit$model$reml,
			moreParams=moreParams, ...)
	
	dimnames(hess) = list(names(baseParam),names(baseParam))
	infmat = solve(hess)*2
	
	pvec = grep("^ci([[:digit:]]|\\.)+$", colnames(fit$summary),
			value=TRUE)
	pvec = as.numeric(gsub("^ci","", pvec))
	qvec = qnorm(pvec)
	names(qvec) = paste("ci", pvec, sep="")
	
	stdErr = sqrt(diag(infmat))
	toAdd = outer(stdErr, qvec, FUN="*")
	
	forSummary = fit$param[colnames(infmat)] + 
			toAdd
	
	
	
	summary = fit$summary
	
	if(any(rownames(forSummary)=='nugget'))
		forSummary = rbind(forSummary,
				sdNugget = sqrt(
						pmax(0,forSummary["nugget",]))
		)
	if(any(rownames(forSummary)=='variance'))
		forSummary = rbind(forSummary,
				sdSpatial = sqrt(
						pmax(0,forSummary["variance",]))
		)
	
	
	inBoth = intersect(rownames(summary), rownames(forSummary))
	
	summary[inBoth,colnames(forSummary)] = 
			forSummary[inBoth,]
	
	list(summary=summary,information=infmat)
}
