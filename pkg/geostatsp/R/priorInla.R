priorInla = function(x, family='gaussian', cellSize=1) {

	names(x) = gsub("^sdNugget$", "sdObs", names(x))

	for(D in names(x)) {
		if(!is.list(x[[D]])) {
			if(length(x[[D]])==1) {
				# pc prior, median is specified
				x[[D]] = c(u=unname(x[[D]]), alpha = 0.5)
			}
			if(is.null(names(x[[D]]))) {
				if(x[[D]][2] > 1) {
					names(x[[D]]) = c('lower', 'upper')
				} else {
					names(x[[D]]) = c('u', 'alpha')
				}
			}
		}
	}

	sdNames = unique(c("sd",grep("^sd", names(x), value=TRUE)))
	if(family=="gaussian") {
		sdNames = unique(c(sdNames, "sdObs"))
	}

    # if model is Gaussian, look for prior for sdNugget
	if(family=="gamma") {
		sdNames = unique(c(sdNames, "gammaShape"))
	}
	if(family %in% c("weibull", "weibullsurv") ) {
		sdNames = unique(c(sdNames, "weibullShape"))
	}

	ciTarget = log(c(0.025, 0.975))
	parGetRid = c('cellSize','info','dprior')


	precPrior=list()

	for(Dsd in sdNames) {
		Dprec = gsub("^sd","precision",Dsd)

		if(is.null(x[[Dsd]])) { 
		# no prior supplied
          # default prior
			x[[Dsd]] = = c(u=1, alpha = 0.5)
		} 
		if(is.list(x[[Dsd]])) { 
	# use it unchanged
			precPrior[[Dsd]] = x[[Dsd]]
			if(!is.null(precPrior[[Dsd]]$initial))
				precPrior[[Dsd]]$initial = precPrior[[Dsd]]$initial^(-2)
		} else if(is.matrix(x[[Dsd]])) {

  # if it's a matrix first column is sd, 2nd column is density

			precPrior[[Dsd]]$dprior = approxfun(
				x[,1], x[,2]
				)

  # get rid of zero
			x[[Dsd]] = x[[Dsd]][x[[Dsd]][,1]>0,]
			logPrecDens = cbind(
				-2 * log(x[[Dsd]][,1]),
				x[[Dsd]][,2] * x[[Dsd]][,1]/2
				)
  # reorder smallest to largest precision		
			logPrecDens = logPrecDens[order(logPrecDens[,1]), ]

			precPrior[[Dsd]] = x[[Dsd]]
			precPrior[[Dsd]]$string = paste0(
				"prior='table:", 
				paste(as.vector(logPrecDens), collapse=' '), 
				"'", sep='')
		} else if(all(c('u','alpha') %in% names(x[[Dsd]]))) {
              # pc priors
			expRate = -log(x[[Dsd]]['alpha'])/x[[Dsd]]['u']*log(2)
			precPrior[[Dsd]] = list(
				param = x[[Dsd]],
				prior = 'pc.prec',
				rprior = eval(parse(text=paste0('function(x) dexp(x, rate=', 1/expRate, ')')))
				)
		} else if(all(c('lower','upper') %in% names(x[[Dsd]]))) {

              # gamma prior for sd
# quantiles provided
			cifun = function(pars) {
				theci = pgamma(
					x[[Dsd]],
					shape=pars[1], 
					rate=pars[2],log.p=TRUE)

				sum(c(1,4)*(ciTarget - theci)^2)
			}

			precPrior2=stats::optim(c(.5,.5/mean(x[[Dsd]])), 
				cifun, 
				lower=c(0.000001,0.0000001),method="L-BFGS-B")


# sd is gamma
# log(sd) is loggamma
# log(prec) = -2*log(sd) is dist'n below

# a=1/3;b=3/4;stuff = rgamma(10000, a,b)
# hist(log(stuff^(-2)), prob=TRUE)
# xseq = seq(-20, 20, len=1001) 
# lines(xseq, 
#  exp(a*log(b)  - lgamma(a) -log(2) - 
#  (a-1)*(xseq/2)-b*exp(-xseq/2)-xseq/2))

			names(precPrior2$par) = c("shape","rate")
			precPrior[[Dsd]]$param = precPrior2$par
			precPrior[[Dsd]]$info = 'gamma prior for sd'
			precPrior[[Dsd]]$string = paste0(
				"prior=\"expression:
				a = ", precPrior[[Dsd]]$param['shape'], ";
				b = 1.0/", precPrior[[Dsd]]$param['rage'], ";
				return (a - 1.0) * (log_precision / 2.0) - 
				b*exp(-log_precision / 2.0) - log_precision/2.0;\"", 
				sep='')
			precPrior[[Dsd]]$extra = list(
				userPriorCI = NULL,
				priorCI = NULL,
				paramUser = precPrior[[Dsd]]$param * 1,
				)
		}

		precPrior[[Dsd]]$rprior = eval(parse(text=paste0('function(x) dgamma(x, shape=',
			precPrior[[Dsd]]['shape'],',rate=', precPrior[[Dsd]]['range'], ')')))

		if(!length(precPrior[[Dsd]]$string)) {
			precPrior[[Dsd]]$string = deparse(
				precPrior[[Dsd]][setdiff(names(precPrior[[Dsd]]), parGetRid], 
					control=NULL)
			}

		} # end loop Dsd

# range
		rangePrior = NULL
		if("range" %in% names(x)) {
			if(is.null(x$range)) 
				x = c(u=10*cellSize, alpha=0.5) 

			if(is.list(x$range)) {
				if(!is.null(x$range$cellSize))
					x$range$cellSize = cellSize
				if(!is.null(x$range$initial))
					x$range$initial = x$range$cellSize/x$range$initial
			}

			if(is.matrix(x$range)) {

				x$range = x$range[x$range[,1]>0,]

          # prior on log(1/(range/cellsize)) = log(cellsize) - log(range)
				logRangeDens = cbind(
					log(cellSize)-log(x$range[,1]),
					x$range[,2] * x$range[,1])

          # reorder smallest to largest
				logRangeDens = logRangeDens[nrow(logRangeDens):1, ]
				rangePrior = list(string=paste0(
					"prior='table: ", 
					paste0(as.vector(logRangeDens), collapse=' '),
					"'", sep=''),
				cellSize = cellSize,
				info = 'table provided for range')
			} else if(
				all(c('u', 'alpha') %in% names(x$range) ) ) { 
    # pc prior
				rangePrior = pcPriorRange(
					q = x$range['u'],
					p = x$range['alpha'],
					cellSize = cellSize
					)
			} else {
		        # gamma prior for scale
				scaleTarget = sort(cellSize/x$range[c('lower','upper')])
				cifun = function(pars) {
					theci = pgamma(
						scaleTarget,
						shape=pars[1], 
						rate=pars[2])

					sum(c(1,4)*(ciTarget - theci)^2)
				}

				scalePrior2=stats::optim(
					c(.5,.5/mean(scaleTarget)), 
					cifun, 
					lower=c(0.000001,0.0000001),method="L-BFGS-B")

				names(scalePrior2$par) = c("shape","rate")
				rangePrior$param = scalePrior2$par
				rangePrior$info = 'gamma prior for scale'

				rangePrior$string = paste0(
					"prior=\"expression:
					a = ", rangePrior$param['shape'], ";
					b = 1.0/", rangePrior$param['rage'], ";
					return (a - 1.0) * (log_range) - 
					b*exp(-log_range) - log_range;\"", 
					sep='')

				rangePrior$dprior = list(range = eval(parse(text=paste0(
					'function(x) x^(-2)*dgamma(1/x, shape=',
					rangePrior$param['shape'],',rate=', rangePrior$param['range'], ')'))),
				scale = eval(parse(text=paste0(
					'function(x) dgamma(x, shape=',
					rangePrior$param['shape'],',rate=', rangePrior$param['range'], ')')))
				)

			}
			if(!length(rangePrior$string))
				rangePrior$string = deparse(rangePrior[setdiff(names(rangePrior), parGetRid), control=NULL)
		}

      # prior for gamma shape
      # log-normal, priorCI is 4 standard deviations

		familyShapeName = grep(
			"(family|gamma|weibull)Shape", names(x), value=TRUE)
		if( length(familyShapeName) ) {
			familyShapeName = familyShapeName[1]
			familyShapePrior  = list(
				prior='gaussian',
				param=c(
					mean=as.numeric(mean(log(x[[familyShapeName]]))),
					precision = as.numeric(abs(diff(log(x[[familyShapeName]])))[1]/4)^(-2)
					)
				)
		} else {
			familyShapePrior = NULL
		}

# to do: function to evaluate priors

		list(
			family = familyShapePrior,
			prec = precPrior,
			rate = ratePrior
			)

	}
