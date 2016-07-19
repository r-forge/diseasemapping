#' @export
precToSd = function(densmat) {
	densmat[,"y"] = densmat[,"y"] * 2*  densmat[,"x"]^(3/2) 
	densmat[,"x"] = 1/sqrt(densmat[,"x"])  
	densmat
}

#' @export
priorPostSd = function(
		res, 
		param=1:length(res$all.hyper$random),
		minSd = 0.001
){
	
	# return prior and posterior of all standard deviation parameters
	
	result = list()
	
	names(res$all.hyper$random) = 
			unlist(lapply(res$all.hyper$random, function(qq) qq$hyperid))
	
	param = names(res$all.hyper$random[param])
	
	thesummary = 1/sqrt(res$summary.hyperpar[paste("Precision for", param),])
	quantCols = grep("quant$", colnames(thesummary))
	Squant = 1-as.numeric(gsub("quant", "", colnames(thesummary)[quantCols]))
	thesummary[,quantCols] = thesummary[,quantCols[order(Squant)]]
	
	for(Dparam in param){
		
		
		Dprec = paste("Precision for", Dparam)
		
		postPrec = res$marginals.hyperpar[[Dprec]]
		
		if(!is.null(postPrec)) {
			result[[Dparam]] = list()
			thesummary[Dprec,'mean'] = INLA::inla.emarginal(
					function(qq) 1/sqrt(qq),
					postPrec
			)
			
			thesummary[Dprec,'sd'] = sqrt(INLA::inla.emarginal(
							function(qq) 1/qq,
							postPrec
					) - thesummary[Dprec,'mean']^2)
			
			result[[Dparam]]$posterior = precToSd(postPrec)
			
			xhyper = res$all.hyper$random[[Dparam]]$hyper
			xhyper = xhyper[[which(unlist(lapply(xhyper, function(qq) qq$short.name)) == 'prec')[1] ]]
			
			
			if(!all(xhyper$prior == 'pc.prec'))
				warning("pc.prec prior expected, got ", xhyper$prior[1])
			
			xSeq = seq(
					minSd, 
					2*max(result[[Dparam]]$posterior[,'x']), 
					len=1000)
			pSeq = xSeq^(-2)
			
			
			result[[Dparam]]$prior = precToSd(
					cbind(
							x=pSeq, 
							y=INLA::inla.pc.dprec(
									pSeq,
									u=xhyper$param[1], 
									alpha=xhyper$param[2]
							)
					) 
			)
		}
		
		# spatial proportion
		Dphi = paste("Phi for", Dparam)
		postPhi = res$marginals.hyperpar[[Dphi]]
		
		if(!is.null(postPhi)) {
			
			priorProp = matrix(scan(text=gsub("[[:alpha:]]+:", "", 
									res$all.hyper$random[[Dparam]]$hyper$theta2$prior),
							quiet=TRUE), 
					ncol=2, dimnames = list(NULL, c( "xTrans","logDensTrans")))
			priorProp = cbind(priorProp, logX = exp(priorProp[,'xTrans']))
			priorProp = cbind(priorProp, 
					x = priorProp[,'logX']/ (1+priorProp[,'logX'])
			)
			
			# what the prior integrates to	
			const = sum(exp(priorProp[,'logDensTrans']+
									c(NA,log(abs(diff(priorProp[,'xTrans']))))),
					na.rm=TRUE)		
			# add to it integrates to 1		
			
			diffTrans = c(NA, diff(priorProp[,'xTrans']))
			diffX = c(NA, diff(priorProp[,'x']))
			
			priorProp = cbind(priorProp, 
					logDens = priorProp[,'logDensTrans'] +
							log(abs(diffTrans)) - 
							log(abs(diffX)) - log(const)
			)
			priorProp[1, 'logDens'] = priorProp[2, 'logDens']
			
			priorProp = cbind(priorProp, 
					dens = exp(priorProp[,'logDens'] )
			)
			priorProp = cbind(priorProp, 
					cDens = cumsum(priorProp[,'dens'] *
  								c(0, abs(diff(priorProp[,'x']))))
			)
			
			result[[Dphi]]$prior = priorProp[,c('x','dens')]
			result[[Dphi]]$posterior = postPhi
		}
	}
	
	if(length(result)==1) result = result[[1]]
	rownames(thesummary) = gsub("^Precision", "SD", rownames(thesummary))
	result$summary = thesummary
	
	result
}


