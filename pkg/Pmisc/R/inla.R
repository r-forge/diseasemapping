#' precisions to standard deviations
#' 
#' @aliases priorPost priorPostSd precToSd
#' @description Transforms prior and posterior distributions of precision parameters to standard deviations
#'
#' @param object an inla result
#' @param param vector of parameters to transform
#' @param group random effect parameters or 'family' parameters
#' @param densmat matrix with columns 'x' and 'y'
#' @export
#' 
priorPostSd = function(
  object, 
  param,
  group = c('random','family')
){
  group = group[1]
  if(missing(param)) {
  param = 1:length(object$all.hyper[[group]])
  }
  param = param[ param <= length(object$all.hyper[[group]]) ]

	 # return prior and posterior of all standard deviation parameters

  res = object
  result = list()
  names(res$all.hyper[[group]]) = 
  unlist(lapply(res$all.hyper[[group]], function(qq) qq$hyperid))

  paramIndex = param
  param = names(res$all.hyper[[group]][paramIndex])

  if(group == 'random') {
    paramLong = paste("Precision for", param)
  } else {
    Slabel = lapply(res$all.hyper[[group]], function(qq) c(
      unlist(lapply(qq$hyper, function(qqq) qqq$short.name)),
      qq$label)
    )
    Slabel = do.call(rbind, Slabel)
    if(ncol(Slabel)==3) Slabel = Slabel[,-2, drop=FALSE]

    Slabel[,2] = gsub("^t$", 'student-t', Slabel[,2])
    paramLong1 = paste0(Slabel[,1], ' (parameter )?for ', Slabel[,2])
    paramLongRegexp = gsub("prec ", "[pP]rec(ision)? ", paramLong1)
    paramLongRegexp = gsub("[?]for ", "?for( the)? ", paramLongRegexp)
    paramLong = unlist(mapply(
      grep, 
      pattern = paramLongRegexp, 
      MoreArgs = list(x=rownames(res$summary.hyper), 
      value=TRUE, ignore.case=TRUE)))
    if(!length(paramLong)) {
      paramLong = grep(paste(
        "^precision for the ", Slabel[,ncol(Slabel)], ' observations$', sep=''
        ), rownames(res$summary.hyper), 
      value=TRUE, ignore.case=TRUE)
    }
  }
  names(paramLong) = param
  
  thesummary = 1/sqrt(res$summary.hyperpar[paramLong,])
  
  quantCols = grep("quant$", colnames(thesummary))
  Squant = 1-as.numeric(gsub("quant", "", colnames(thesummary)[quantCols]))
  thesummary[,quantCols] = thesummary[,quantCols[order(Squant)]]
  thesummary[,'mode'] = NA

  Slty = c(prior=2, posterior=1)
  Scol = c(prior='red', posterior='black')
  
  for(Dparam in param){
    Dprec = paramLong[Dparam]
 
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

      thesummary[Dprec,'mode'] = result[[Dparam]]$posterior[
        which.max(result[[Dparam]]$posterior[,2]), 1]

      xhyper = res$all.hyper[[group]][[
      which(Dparam==param)]]$hyper[[1]]
#			   xhyper = xhyper[[which(unlist(lapply(xhyper, function(qq) qq$short.name)) == 'prec')[1] ]]



      forXlim = 20*c(0.8, 1.2)*range(thesummary[
        Dprec, 
        grep("quant$", colnames(thesummary))])

      forXlim = range(c(floor(forXlim), ceiling(forXlim)))/20

      newXseq = sort(c(forXlim, result[[Dparam]]$posterior[,1]))

      maxX = 2*min(res$summary.hyperpar[Dprec,
        grep('quant$', colnames(res$summary.hyperpar))])^(-1/2)


      xSeq = seq(
        0, 
        maxX, 
        len=1000)			   

      result[[Dparam]]$posterior = cbind(
        x = newXseq,
        y = stats::approx(
          result[[Dparam]]$posterior[,1],
          result[[Dparam]]$posterior[,2],
          newXseq, rule=2
          )$y)


      if(all(trimws(xhyper$prior) == 'loggamma')) {

        result[[Dparam]]$prior = 
        cbind(
          x=xSeq, 
          y=stats::dgamma(
            xSeq,
            rate = xhyper$param[2],
            shape = xhyper$param[1]
            )
          )

        result[[Dparam]]$posterior = cbind(
          result[[Dparam]]$posterior,
          prior = stats::dgamma(
            newXseq,
            rate = xhyper$param[2],
            shape = xhyper$param[1]
            )
          )
      } else if(all(trimws(xhyper$prior) == 'pc.prec')) {

        result[[Dparam]]$prior = 
        cbind(
          x=xSeq, 
          y=stats::dexp(
            xSeq,
            -log(xhyper$param[2])/xhyper$param[1]
            )
          )

        result[[Dparam]]$posterior = cbind(
          result[[Dparam]]$posterior,
          prior = stats::dexp(
            newXseq, 
            -log(xhyper$param[2])/xhyper$param[1]
            )
          )

      } else {
        warning("pc.prec prior expected, got ", xhyper$prior[1])
      } # end if statements

      if( (forXlim[1]<0.25) & (forXlim[2]>0.5)) forXlim[1] = 0

      result[[Dparam]]$matplot = list(
        x = result[[Dparam]]$posterior[,1], 
        y = result[[Dparam]]$posterior[,3:2], 
        type='l', lty=Slty, 
        col=Scol, 
        xlab='sd',
        ylab='dens', 
        xaxs='i',
        xlim = forXlim      )

    } # postPrec not null

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

      xLim = rev(c(0.8, 1.2))*range(result[[Dphi]]$posterior[,1])
      if( (xLim[1]<0.25) & (xLim[2]>0.5)) xLim[1] = 0

      result[[Dphi]]$matplot = list(
        x = result[[Dphi]]$posterior[,1], 
        y = result[[Dphi]]$posterior[,3:2], 
        type='l', lty=Slty, 
        col=Scol, 
        xlab='sd',
        ylab='dens', 
        xaxs='i',
        xlim = xLim
        )


    } # not null postphi
  } # loop through params

  if(length(result)==1) result = result[[1]]
  rownames(thesummary) = gsub("^[Pp]recision", "SD", rownames(thesummary))
  result$summary = thesummary

  result$legend = list(
    x="topright", lty=Slty, 
    col=Scol, 
    legend=names(Slty), 
    bty='n')

  result
}


#' @export
precToSd = function(densmat) {
  densmat[,"y"] = densmat[,"y"] * 2*  densmat[,"x"]^(3/2) 
  densmat[,"x"] = 1/sqrt(densmat[,"x"])  
  densmat
}


#' @export
priorPost = function(object) {

  Slty = c(prior=2, posterior=1)
  Scol = c(prior='red', posterior='black')

  paramList = lapply(
    object$all.hyper[c('family','random')],
    function(x) {
      resx = lapply(x, # apply over parameters withing group
        function(xx) {
          resxx = lapply(xx$hyper, # apply over hyperparameters for this prior
            function(xxx) {
              unlist(xxx[c('name','short.name','prior')])
            })
          resxx = do.call(cbind, resxx)
          resxx = rbind(resxx, 
            internal.name = colnames(resxx)
            )
          colnames(resxx) = resxx['short.name',]
          resxx = rbind(resxx, id=xx$hyperid,
            label = c(xx$label, NA)[1])
          resxx
        })
      Sindex = rep(1:length(resx), unlist(lapply(resx, ncol)))
      resx = do.call(cbind, resx)
      if(is.null(colnames(resx))) colnames(resx) = resx['label', ] 
      colnames(resx) = paste(colnames(resx), resx['id',], sep='_')
      resx = as.data.frame(t(resx))
      resx$index = Sindex
      resx
    })

  notInFamily = setdiff(names(paramList$random), names(paramList$family))
  paramList$family = cbind(
    paramList$family,
    matrix(NA, 
      nrow = nrow(paramList$family),
      ncol = length(notInFamily),
      dimnames = list(rownames(paramList$family), notInFamily)
      )
    )[, colnames(paramList$random)]
  paramDf = do.call(rbind, paramList)
  paramDf$group = rep(names(paramList), unlist(lapply(paramList, nrow)))
  paramDf$out.name = NA
  paramDf = paramDf[grep('precision offset', paramDf[,'name'], invert=TRUE), , drop=FALSE]


  # a few fixes
  # weibull cure has short name a which is alpha in summary table
  paramDf$short.name = gsub("^a$", "alpha", paramDf$short.name)
  paramDf$prior = trimws(paramDf$prior)

  result = list()
  for(Dparam in which(!is.na(paramDf$prior))) {
    if(paramDf[Dparam, 'prior'] == 'pc.prec') {
      result[[Dparam]] = priorPostSd(
        object, 
        param=paramDf[Dparam, 'index'],
        group = paramDf[Dparam, 'group']
        )
      paramDf[Dparam, 'out.name'] = paste("sd for", 
        gsub("INLA.Data", paramDf[Dparam, 'label'], paramDf[Dparam, 'id']))
    } else {
      # not a pc prec prior
      paramDf[Dparam, 'out.name'] = paste(
        paramDf[Dparam, 'short.name'], 'for',
        gsub("INLA.Data", paramDf[Dparam, 'label'], paramDf[Dparam, 'id']))

      if(paramDf[Dparam, 'short.name'] == 'prec') {
        inlaNameHere = paste(
          'Precision for', paramDf[Dparam, 'id']
          )  
      } else {
        paramName = gsub("^weibullcure$", "ps", 
          as.character(stats::na.omit(unlist(paramDf[Dparam, c('label','id')][1])))
          )
        inlaNameHere = 
        grep(paste(
          '^', 
          gsub("^prob$", "p", paramDf[Dparam, 'short.name']), 
          '( parameter)? for ', paramName, 
          '$',
          sep=''),
        rownames(object$summary.hyperpar),
        ignore.case=TRUE, value=TRUE)
      }
      # a fix for weibullcure models
      inlaNameHere = gsub("weibullcure", "ps", inlaNameHere)
      inlaNameHere = gsub("^prob parameter for ps", "p parameter for ps", inlaNameHere)

      priorHere = object$all.hyper[[
      paramDf[Dparam, 'group'] 
      ]][[
      paramDf[Dparam, 'index']  
      ]][['hyper']][[  
      paramDf[Dparam, 'internal.name']  
      ]]


      postHere = object$marginals.hyperpar[[inlaNameHere]]

      xLim = range(object$summary.hyperpar[inlaNameHere,
        grep('quant$', colnames(object$summary.hyperpar))
        ])*c(0.7, 1.5)*20
      xLim = range(c(floor(xLim), ceiling(xLim)))/20
      if( (xLim[1]<0.25) & (xLim[2]> (2*xLim[1]) ) ) xLim[1] = 0

      if(is.null(postHere)) {
        xRange = xLim
        xSeq = seq(xRange[1], xRange[2], len=1000)
        postHere = cbind(
          x = xSeq, y = NA
          )
      } else {
        xSeq = postHere[,1]
        xRange = range(xSeq)
      }

      if(priorHere$prior=='normal'){
        priorDens = stats::dlnorm(
          xSeq, meanlog = priorHere$param[1],
          sdlog = priorHere$param[2]
          )

        xRangeP = stats::qlnorm(
          c(0.01, 0.99), 
          meanlog = priorHere$param[1],
          sdlog = priorHere$param[2]
          )
        xSeqP = seq(xRangeP[1], xRangeP[2], len=1000)
        priorMat = cbind(x=xSeqP, 
          y=stats::dlnorm(
            xSeqP, meanlog = priorHere$param[1],
            sdlog = priorHere$param[2]
            ))

      } else if(priorHere$prior == 'loggamma') {

        priorDens = stats::dgamma(xSeq, 
          shape=priorHere$param[1], 
          rate=priorHere$param[2])

        xRangeP = stats::qgamma(
          c(0.001, 0.999), 
          shape=priorHere$param[1], 
          rate=priorHere$param[2])

        xSeqP = seq(xRangeP[1], xRangeP[2], len=1000)

        priorMat = cbind(x=xSeqP, 
          y=stats::dgamma(
            xSeqP,           
            shape=priorHere$param[1], 
            rate=priorHere$param[2])
          )

      } else if(priorHere$prior == 'logitbeta') {
        priorDens = stats::dbeta(xSeq, 
          shape1=priorHere$param[1], 
          shape2=priorHere$param[2])

        xRangeP = stats::qbeta(
          c(0.001, 0.999), 
          shape1=priorHere$param[1], 
          shape2=priorHere$param[2])
        xSeqP = seq(xRangeP[1], xRangeP[2], len=1000)
        priorMat = cbind(x=xSeqP, 
          y=stats::dbeta(xSeqP, 
            shape1=priorHere$param[1], 
            shape2=priorHere$param[2])
          )

      } else {
        warning(priorHere$prior, 'not implemented yet')
        priorDens = NA
        priorMat = NULL
      }

      result[[Dparam]] = list(
        prior = priorMat,
        posterior = cbind(postHere, prior=priorDens)
        )

      result[[Dparam]]$matplot = list(
        x = result[[Dparam]]$posterior[,1], 
        y = result[[Dparam]]$posterior[,3:2], 
        type='l', lty=Slty, 
        col=Scol, 
        xlab=gsub(" .*$", "", paramDf[Dparam, 'out.name']),
        ylab='dens', 
        xaxs='i',
        xlim = xLim
        )

      result[[Dparam]]$summary = object$summary.hyperpar[inlaNameHere,]      
    } # parameter not pc.prec
  } # loop through parameters

  names(result) = paramDf[,'out.name']

  result$parameters = paramDf[,'out.name']

  result$legend = list(
    x="topright", lty=Slty, 
    col=Scol, 
    legend=names(Slty), 
    bty='n')

  sdSeq = grep("^sd ", names(result), value=TRUE)
  if(length(sdSeq)) {
    theLims = lapply(result[sdSeq], function(x) cbind(
      range(x$matplot$xlim),
      range(x$matplot$y)))
    theLims = do.call(rbind, theLims)
    theLims = apply(theLims, 2, range)
    colnames(theLims) = c('x', 'y')
    rownames(theLims) = c('lower','upper')
    theLims['lower','y'] = 0
    result$extra$limits = theLims
    for(D in sdSeq) {
      if(sum(abs(result[[D]]$matplot$xlim - theLims[,'x'])) < 0.1) {
        result[[D]]$matplot$xlim = theLims[,1]
        result[[D]]$matplot$ylim = theLims[,2]
      }
    }
  }

  result$summary = do.call(
    rbind,
    lapply(result[result$parameters], function(x) x$summary)
    )

  if(!length(grep("2$", rownames(result$summary)))) {
    rownames(result$summary) = gsub("1$", "", rownames(result$summary))
  }

  result
}
