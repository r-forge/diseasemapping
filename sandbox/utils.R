#+ precToSdFunction, include=FALSE
precToSd = function(densmat) {
	densmat[,"y"] = densmat[,"y"] * 2*  densmat[,"x"]^(3/2) 
	densmat[,"x"] = 1/sqrt(densmat[,"x"])  
	densmat
}

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

#'


# make multiple plots as pseudo-subfigures, in a table
    hook_plot_p = function(x, options) {
    
    # for debugging
    # stuff <<- list(x=x, options=options)
    
    if (options$fig.show == 'animate') return(hook_plot_html(x, options))
    
    base = opts_knit$get('base.url') 
    if(is.null(base)) base=''
    cap = options$fig.cap
    scap = options$fig.subcap
    if(is.null(cap)) cap=''
    if(is.null(scap)){
    scap = cap
    }
    
    if (is.null(w <- options$out.width) & is.null(h <- options$out.height) &
    is.null(s <- options$out.extra) & options$fig.align == 'default') {
    result = 
    sprintf('![%s](%s%s) ', scap, base, knitr:::.upload.url(x))
    
    } else {
    # use HTML syntax <img src=...>
    result = 
    knitr:::.img.tag(knitr:::.upload.url(x), w, h, scap,
    c(s, sprintf('style="%s"', knitr:::css_align(options$fig.align)))
    )
    }
    
    if(any(options$fig.ncol==0)){
      return(result)
    }
    
    
    fig.num =options$fig.num
    if(is.null(fig.num))
    fig.num=1
    
    
    fig.subcap.all = options$params.src
    fig.subcap.all = eval(parse(text=
    paste("list(", sub("^[[:alnum:]]+,", "", fig.subcap.all), ")")
    ))$fig.subcap
    if(length(fig.subcap.all) < fig.num)
    fig.subcap.all = c(
    fig.subcap.all, 
    rep(" ", fig.num - length(fig.subcap.all))
    )
    fig.subcap.all = fig.subcap.all[1:fig.num]
    
    
    if(length(fig.subcap.all)) {
    fig.cur =options$fig.cur
    if(is.null(fig.cur))
    fig.cur=1
    fig.ncol = options$fig.ncol
    if(is.null(fig.ncol)) {
    fig.ncol = 1
    }
    
    Drow = floor((fig.cur-1)/fig.ncol)+1
    Dcol = fig.cur - (Drow-1) * fig.ncol 
    
#	cat("\n rc ", fig.cur, Drow, " ", Dcol, "\n")
    
    if(Dcol==1) {
    result = paste("|", result, sep="")
    }
    
    if( (Dcol == fig.ncol ) | (fig.cur == fig.num) ) { 
    # we're at the end of a column or the last plot
    Dend = fig.ncol*Drow
    nextra = Dend - fig.cur
    # pandoc doesn't like single column tables, add an extra column
    if(fig.ncol==1) nextra = nextra + 1
    result = paste( result,
    paste(rep("   |  ", nextra), collapse=""),		
    "|\n")
    # add the subcaptions	
    
    Ssubcap = seq(fig.ncol*(Drow-1)+1,
    min(c(Dend, length(fig.subcap.all)))
    )
    fig.subcap.all = fig.subcap.all[Ssubcap]
    fig.subcap.all = paste(letters[Ssubcap], fig.subcap.all, sep=") ")
    fig.subcap.all= c(	fig.subcap.all,
    rep(" ",nextra)
    ) 
    
    result = paste(result, "| ",
    paste(fig.subcap.all, 
    collapse=" | "), 
    " |\n", sep=""
    )
    
    } else {
    result = paste( result, "|", sep="")  	
    }
    
    if(fig.cur==1) {
    result = paste('\n\n\n|',
    paste(rep('     ', max(c(2,fig.ncol))), collapse='|'), '|\n', 
    "|",
    paste(rep(':---:', max(c(2,fig.ncol))), collapse='|'),'|\n', 
    result, sep="")
    }
    
    if(fig.cur == fig.num)
    result = paste(result, "\n : Figure: ", cap, "\n\n\n\n", sep="")	
    
    }
    result
    
    }
    