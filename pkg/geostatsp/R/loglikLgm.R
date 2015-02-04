loglikLgm = function(param,
		data, formula, coordinates=data,
		reml=TRUE, 
		minustwotimes=TRUE,
		stored=NULL, moreParams=NULL) {

	# create 'covariates', and 'observations'
 	
	trend = formula
	
	if(class(trend)=="formula") {
		observations = all.vars(trend)[1]
		if(!any(names(data)==observations))
			warning("can't find observations ", observations, "in data")
		# frame first, so we see which row names have been omitted due to NA's
		# sometimes model.matrix removes row names
		covariates = model.frame(trend, data.frame(data))
		observations = covariates[,	observations]
		theRowNames = rownames(covariates)
		covariates = model.matrix(trend,covariates)
		rownames(covariates) = theRowNames
		theNA = which(!rownames(data.frame(data)) %in% theRowNames)
	} else if(!is.null(trend)){
		# observations must be a vector and covariates a matrix
		covariates= as.matrix(trend)
		theNA = apply(covariates, 1, function(qq) any(is.na(qq)))
		observations=data[!theNA]
		theNA = which(theNA)
	} else {
		theNA = NULL
	}
	
		if(length(grep("SpatialPoints", class(coordinates)))) {
			if(length(theNA))
				coordinates = coordinates[-theNA,]
			
		} else if(	class(coordinates) == "dist")	{
			if(length(theNA)) {
				coordinates = 
					as.matrix(coordinates)[-theNA,-theNA]
      } else {
        coordinates= as.matrix(coordinates)
      } 
    } else if(class(coordinates) == "matrix")	{
      if(length(theNA)) {
        coordinates = coordinates[-theNA,-theNA]
      }
    } else {
				warning("coordinates must be a SpatialPoints object\n or matrix a dist object.  It's being assumed \n it's a matrix of coordinates")
				coordinates = SpatialPoints(coordinates)
				if(length(theNA))				
					coordinates = coordinates[-theNA,]
		}

		# sort out the parameters
		param = c(param, moreParams)
		
		names(param) = gsub("^var$", "variance", names(param))
		
		param=param[!is.na(param)]
	
		haveVariance = any(names(param)=="variance")
		haveNugget = any(names(param)=="nugget")
		if(haveNugget){
			if(param["nugget"] == 0) {
				haveNugget=FALSE
			}
		} else {
			param["nugget"] = 0
		}

		if(!haveVariance) {
			# variance will be estimated.  
		  # nugget is assumed to be nugget/variance
			param["variance"] = 1
		}
			
		# box cox transform
		twoLogJacobian = 0
		if(any(names(param)=="boxcox")) {
			
			if(any(names(stored)=="boxcox")) {
				if(abs(stored$boxcox-param["boxcox"])>0.0001)
					warning("boxcox param and stored are ", param["boxcox"],
							" and ", stored$boxcox)
				twoLogJacobian = stored$twoLogJacobian
				observations = stored$observations
			} else { # need to compute boxcox stuff
			
			if(abs(param["boxcox"]-  1 ) < 0.001) {
			 # boxcox close to 1, don't transform
				twoLogJacobian=0
				
			} else { # box cox is not one.
				twoLogJacobian = 2*(param["boxcox"]-1)* 
					sum(log(observations))	
			
				if(is.nan(twoLogJacobian))
					warning("boxcox shouldnt be used with negative data")

				if(abs(param["boxcox"])<0.001) {
					observations = log(observations) 
				} else  { #boxcox far from 0 and 1
				observations <- 
						((observations^param["boxcox"]) - 1)/
							param["boxcox"]
				}
			} # end boxcox param far from 1
		} # end need to compute boxcox
		} # end have box cox
		
		
    obsCov = cbind(observations, covariates)
    if(FALSE){ # old code
    fromC = maternCholSolve(param, obsCov, coordinates)

    betaHat = fromC$betaHat
    
    useMl = c('ml', 'reml')[1+reml]
    useVariance = c('estimatedVariance','fixedVariance')[1+haveVariance]
    minusTwoLogLik = 
        (fromC$minusTwoLogLik[useMl,useVariance] - twoLogJacobian)
    varBetaHat = fromC$varBetaHat[,,useVariance, useMl]
    totalVarHat = fromC$totalVarHat[useMl]
    
    if(haveVariance) {
      varMle = fromC$varMle[useMl,]
    } else {
      varMle = param[colnames(fromC$varMle)]
      param[c("variance","nugget")] = 
          totalVarHat * param[c("variance","nugget")]
    }
    
  # format the output
  result = minusTwoLogLik
	
} else {
  
  Nobs = nrow(obsCov)
  Ncov = ncol(obsCov)-1
  Nrep = 1
  
  paramFull = fillParam(param)
  Ltype = c(ml=0, reml=1, mlFixed=2, remlFixed=3)
  Ltype = reml + 2*haveVariance

  if(class(coordinates)=='matrix'){
    xcoord = as.vector(coordinates)
    ycoord = -99
    aniso=FALSE
  } else if(length(grep("^SpatialPoints", class(coordinates)))){
    xcoord=coordinates@coords[,1] 
    ycoord=coordinates@coords[,2]
    aniso=TRUE
  } else {
    warning('coordinates should be SpatialPoints or matrix')
    xcoord=ycoord=aniso=NULL
  }
  
  resultC = .C("maternLogL",
      xcoord=as.double(xcoord), 
      ycoord=as.double(ycoord),
      param=as.double(paramFull[
              c('nugget','variance','range',
                  'shape','anisoRatio', 'anisoAngleRadians')]),
      aniso=as.integer(aniso),
      obsCov = as.double(obsCov),
      N= as.integer(c(Nobs,Nrep,Ncov)),
      boxcox=as.double(-9.9),
      boxcoxType=as.integer(0),
      logL=as.double(rep(-9.9, Nrep+1)),
      totalVarHat=as.double(rep(-9.9, Nrep)),
      betaHat = as.double(rep(-9.9, Ncov)), 
      varBetaHat = as.double(rep(-9.9, Ncov* Ncov)),
      Ltype=as.integer(Ltype)
  )

  totalVarHat = resultC$totalVarHat
  betaHat = resultC$betaHat
  names(betaHat) = colnames(obsCov)[-1]
  
  varBetaHat = new("dsyMatrix", 
      Dim = as.integer(c(Ncov, Ncov)), 
      uplo="L",
      x=resultC$varBetaHat)
  dimnames(varBetaHat) = list(names(betaHat),names(betaHat))
  varBetaHat = varBetaHat*totalVarHat

  
  result = resultC$logL[Nrep+1]- twoLogJacobian
      
}

param[c("variance","nugget")] = 
    totalVarHat * param[c("variance","nugget")]



if(minustwotimes) {
  names(result) = "minusTwoLogLik"
} else {
  result = -0.5*result
  
  names(result)="logLik"
} 

if(reml)
  names(result) = gsub("Lik$", "RestrictedLik", names(result))

attributes(result)$param = param
attributes(result)$totalVarHat = totalVarHat
attributes(result)$betaHat = betaHat
attributes(result)$varBetaHat = varBetaHat
attributes(result)$reml=reml

  
result
}

 


likfitLgm = function(
		formula, data,
		coordinates=data,
		param=c(range=1,nugget=0,shape=1),
		upper=NULL,lower=NULL, parscale=NULL,
		paramToEstimate = c("range","nugget"),
		reml=TRUE) {

	trend = formula
	
	# for some reason thing break if I remove this next line...
#	stuff = (class(coordinates))
	theproj = proj4string(coordinates)
	
	# check for the variance parameter
	estimateVariance = TRUE
	if(any(paramToEstimate=="variance")) {
		# remove varinace, it's estimated by profile likelihood
		paramToEstimate = paramToEstimate[paramToEstimate != "variance"]
		param = param[names(param)!="variance"]
	} else {
		if(any(names(param)=="variance")){
			estimateVariance = FALSE
#			warning("variance is fixed and not estimated. If this isn't what you wanted remove variance from param")
		}		
	}
	 
	# limits
	lowerDefaults = c(nugget=0,range=0,anisoRatio=0.0001,
			anisoAngleRadians=-pi/2,anisoAngleDegrees=-90,
			shape=0.01,boxcox=-3,variance=0)
	
	upperDefaults= c(nugget=Inf,range=Inf,anisoRatio=Inf,
			anisoAngleRadians=pi/2,anisoAngleDegrees=90,
			shape=100,boxcox=3,variance=Inf)
	

	lowerDefaults[names(lower)]=lower
	
	upperDefaults[names(upper)] = upper
	
	# par scale
	parscaleDefaults = c(range=NA, #range is delt with below
			nugget=1,
			boxcox=1,
			anisoAngleDegrees=60,
			anisoAngleRadians=1,
			anisoRatio=1,
			variance=1,
			shape=1)

#	if(length(grep("^SpatialPoints", class(coordinates))))
#		print("wah")
 	
	# convert input data to a model matrix
	if(class(trend)=="formula") {
 
		data = data.frame(data)
		theNA = apply(
				data[,all.vars(trend)[-1],drop=FALSE],
				1, function(qq) any(is.na(qq)))
		noNA = !theNA
		
		covariates = model.matrix(trend, data[noNA,])
		observations = all.vars(trend)[1]
		
		if(!any(names(data)==observations))
			warning("can't find observations ", observations, "in data")
		observations = data[noNA,observations]
 
		
	} else {
		# observations must be a vector and covariates a matrix
		trend = as.matrix(trend)
		theNA = is.na(data) | apply(trend, 1, 
					function(qq) any(is.na(qq))
				)
		noNA = !theNA
				
		observations=data[noNA]
		covariates=trend[noNA,,drop=FALSE]

	}

	if(any(theNA)) {
		if(length(grep("^SpatialPoints", class(coordinates)))) {
			theRowNames= rownames(data.frame(coordinates))[noNA]
			coordinates = SpatialPoints(coordinates)[noNA]	
		} else if(class(coordinates)=="dist"){
			coordinates = as.matrix(coordinates)
			coordinates = coordinates[noNA,noNA]
			theRowNames = rownames(coordinates)
			coordinates = as.dist(coordinates)
			
		} else {
			theRowNames = NULL
			warning("missing vlaues in data but unclear how to remove them from coordinates")
		}
	} else {
		theRowNames = rownames(data.frame(coordinates))
	}
	
	# if the model's isotropic, calculate distance matrix

 	coordinatesOrig = coordinates
	  # not estimating aniso params and coordinates is SpatialPoints
	if( (!length(grep("^aniso", paramToEstimate))) &
			length(grep("^SpatialPoints", class(coordinates)))) {

		# see if there's anisotropy which is fixed, not estimated
	# which would be odd but we'll test nonetheless
		if(length(grep("^anisoRatio", names(param)))){
			if(abs(param["anisoRatio"]- 1) > 0.0001){
				# it is indeed anisotropic
			
				if(any(names(param)=="anisoAngleDegrees") & 
						!any(names(param)=="anisoAngleRadians") ) {
					param["anisoAngleRadians"] = param["anisoAngleDegrees"]*2*pi/360				
				}
	
				if(class(coordinates)=='dist')
					warning("distance matrix supplied but anisotropic model is used.  Supply points instead.")
				x = coordinates@coords[,1] + 1i*coordinates@coords[,2]
				
				x = x * exp(1i*param["anisoAngleRadians"])
				x = Re(x) +  (1i/ param["anisoRatio"] )*Im(x)
				coordinates = SpatialPoints(cbind(Re(x), Im(x)))
			} # end is anisotripic		
		} # end anisotropy params supplied
		# need to assign names manually
		theRowNames = rownames(data.frame(coordinates))
		coordinates = dist(coordinates(coordinates))
		attributes(coordinates)$Labels = theRowNames
 
	}  
	parscaleDefaults["range"] = dist(t(bbox(coordinatesOrig)))/2
	parscaleDefaults[names(parscale)] = parscale
	
	
	# default starting values for parameters
	paramDefaults = c(nugget=0,anisoRatio=1, anisoAngleDegrees=0,
			anisoAngleRadians=0,shape=1, boxcox=1,
			range=as.numeric(parscaleDefaults["range"]/2))

	
	startingParam = param[paramToEstimate]
	names(startingParam) = paramToEstimate # fixes names lost when no starting value provided

	naStarting = is.na(startingParam)
	startingParam[naStarting]= paramDefaults[names(startingParam)[naStarting]]
	
	moreParams = param[!names(param) %in% paramToEstimate]
	

	# check to see if it's worth storing box cox quantities
if(any(paramToEstimate=="boxcox")) {
	if(any(observations<=0)){
		warning("box cox transform specified with negative observations")
	}
}

if(any(names(param)=="boxcox") ) {
	
	if(param["boxcox"]!= 1)  {
		
		if(any(observations<=0)){
			warning("box cox transform specified with negative observations")
		}
		
	}		
	
	
}		
		
if(any(names(param)=="boxcox") & !any(paramToEstimate=="boxcox")) {

		if(abs(param["boxcox"]-1)>0.0001){ # boxcox not 1
			stored = list(
					boxcox = param["boxcox"],
					twoLogJacobian = 2*(param["boxcox"]-1)* 
							sum(log(observations))	
			)
			if(abs(param["boxcox"])<0.001) {
				stored$observations = log(observations) 
			} else  { #boxcox far from 0 and 1
				stored$observations <- 
						((observations^param["boxcox"]) - 1)/
						param["boxcox"]
			}
			
		} else { # boxcox is 1
			stored=NULL
		}
	} else { # no box cox
		stored=NULL
	}
	
  if(FALSE){
    
  allParams = c(startingParam, moreParams)
  allParams = geostatsp:::fillParam(allParams)
  paramsForC = allParams[c('nugget','variance','range','shape','anisoRatio','anisoAngleRadians','boxcox')]

  Sparam = names(paramsForC) %in% gsub("AngleDegrees$","AngleRadians", paramToEstimate)
  names(Sparam) = names(paramsForC)
  Sparam['boxcox'] = FALSE

  
  if(class(coordinates)=='matrix'){
    xcoord = as.vector(coordinates)
    ycoord = -99
    aniso=FALSE
  } else if(length(grep("^SpatialPoints", class(coordinates)))){
    xcoord=coordinates@coords[,1] 
    ycoord=coordinates@coords[,2]
    aniso=TRUE
  }
  
  obsCov = cbind(y1=observations, y2=0, y3=0, covariates)

  parToEstC = gsub("AngleDegrees$",
      "AngleRadians", paramToEstimate)
  parToEstC = names(Sparam)[names(Sparam) %in% parToEstC]
  
  forO = list(
      scalarF = c(fnscale=-1, 
          abstol=-1,
          reltol = -1,
          alpha=-1, beta=-1,gamma=-1,
          factr= 1e7,
          pgtol=0.1),
      scalarInt=c(trace=0,
          maxit=200,
          REPORT=1,
          type=-1,lmm=12,
          tmax=-1,temp=-1
      ),
      pars = cbind(
          lower=geostatsp:::fillParam(
              lowerDefaults)[parToEstC], 
          upper=geostatsp:::fillParam(
              upperDefaults)[parToEstC],
          parscale = geostatsp:::fillParam(
              parscaleDefaults)[parToEstC],
          ndeps=rep(0.001, length(parToEstC))
    )
  )
  
  forO$parsInt = rep(0, nrow(forO$pars))
  names(forO$parsInt) = rownames(forO$pars)
  forO$parsInt[
      forO$pars[,'lower'] != -Inf &
          forO$pars[,'upper'] != Inf
      ] = 2

  forO$parsInt[
          forO$pars[,'lower'] != -Inf & 
              forO$pars[,'upper'] == Inf
  ] = 1

  forO$parsInt[
      forO$pars[,'lower'] == -Inf & 
          forO$pars[,'upper'] != Inf
  ] = 3
  
  forO$pars[!is.finite(forO$pars)]=-1
  forO$pars = c(forO$pars, rep(0.0, nrow(forO$pars)^2))
  
  fromOptim = .C(
          "maternLogLOpt",
      start=paramsForC,    
      Sparam=as.integer(Sparam),
      obsCov=obsCov, 
      xcoord,
      ycoord,
      as.integer(aniso),
      as.integer(c(nrow(obsCov), 3, ncol(covariates))),
    Ltype=as.integer(reml),
    optInt = as.integer(forO$scalarInt),
    optF = as.double(forO$scalarF),
    betas=cbind(forO$pars),
    limType = as.integer(forO$parsInt),
    message=format(" ",width=80)
  )

  fromOptim = list(
    mle=fromOptim$start,
    betaHat = fromOptim$betas[1:ncol(covariates)],
    m2logL = fromOptim$optF[1],
    totalVarHat = fromOptim$optF[2],
    message = fromOptim$message,
    detail = fromOptim$optInt[1:3],
    varBetaHat =  
        new("dsyMatrix", 
            Dim = as.integer(rep(ncol(covariates),2)), 
            uplo="L",
            x=fromOptim$betas[
                seq(1+ncol(covariates), len=ncol(covariates)^2)]
        )
  )
  
  fromOptim$totalVarHat = 
      fromOptim$totalVarHat * fromOptim$totalVarHat
  names(fromOptim$detail)  = c(
      "fail","fncount","grcount"
        )
 names(fromOptim$betaHat) = colnames(covariates)
 dimnames(fromOptim$varBetaHat) = list(
     names(fromOptim$betaHat),names(fromOptim$betaHat)
     )
     
  fromOptim$mle[c('nugget', 'variance')] = 
      fromOptim$mle[c('nugget', 'variance')] * fromOptim$totalVarHat  
  
p1=.C("maternLogL",
      xcoord=as.double(xcoord), 
      ycoord=as.double(ycoord),
      param=paramsForC,
      aniso=as.integer(aniso),
      obsCov = as.double(obsCov),
      N= as.integer(c(nrow(obsCov), 3, ncol(covariates))),
      boxcox=as.double(rep(-9.9,9)),
      boxcoxType=as.integer(0),
      logL=as.double(rep(-9.9, 3+1)),
      totalVarHat=as.double(rep(-9.9, 3+1)),
      betaHat = as.double(rep(-9.9, ncol(covariates)*3) ), 
      varBetaHat = as.double(rep(-9.9, ncol(covariates)* ncol(covariates))),
      Ltype=as.integer(reml)
  )
  
  p2 = paramsForC
  p2['nugget']=0.01
  p2=.C("maternLogL",
      xcoord=as.double(xcoord), 
      ycoord=as.double(ycoord),
      param=p2,
      aniso=as.integer(aniso),
      obsCov = as.double(obsCov),
      N= as.integer(c(nrow(obsCov), 3, ncol(covariates))),
      boxcox=as.double(rep(-9.9,9)),
      boxcoxType=as.integer(0),
      logL=as.double(rep(-9.9, 3+1)),
      totalVarHat=as.double(rep(-9.9, 3+1)),
      betaHat = as.double(rep(-9.9, ncol(covariates)*3)), 
      varBetaHat = as.double(rep(-9.9, ncol(covariates)* ncol(covariates))),
      Ltype=as.integer(reml)
  )
  cbind(p1$logL,p2$logL)
  apply(cbind(p1$logL,p2$logL),1,diff)
  
  list(temp2$betaHat, temp$betaHat)
  list(temp2$logL, temp$logL)

  

      
  
  tempo=o2(fn=loglikLgm, par=startingParam, 
			lower=lowerDefaults[paramToEstimate], 
			upper=upperDefaults[paramToEstimate],
			control = list(parscale=parscaleDefaults[paramToEstimate]),
			data=observations, 
			formula=covariates, coordinates=coordinates,
		reml=reml, moreParams=moreParams,
		method = "L-BFGS-B")
  
  }
  
	fromOptim = optim(fn=loglikLgm, par=startingParam, 
			lower=lowerDefaults[paramToEstimate], 
			upper=upperDefaults[paramToEstimate],
			control = list(parscale=parscaleDefaults[paramToEstimate]),
			data=observations, 
			formula=covariates, coordinates=coordinates,
		reml=reml, moreParams=moreParams,
		method = "L-BFGS-B", stored=stored
		)
    
    
	fromOptim$start = startingParam
	fromOptim$parscale = parscaleDefaults[paramToEstimate]

		
		
	fromLogLik = loglikLgm(param=fromOptim$par,
			moreParams=moreParams, 
			data=observations, 
			formula=covariates, coordinates=coordinates,
			reml=reml
			)
			
	result = list(
			param=c(attributes(fromLogLik)$betaHat,
					attributes(fromLogLik)$param),
			varBetaHat = attributes(fromLogLik)$varBetaHat,
			optim=fromOptim, 
			data = data.frame(
					observations = observations,
					resid=observations -
              covariates %*% attributes(fromLogLik)$betaHat
			)
	)
	rownames(result$data) = theRowNames
	
	if(any(names(result$param)=="boxcox") ) {
		
		if(abs(result$param["boxcox"]-1)>0.0001){ # boxcox not 1

			if(abs(result$param["boxcox"])<0.001) {
				result$data$obsBC = log(observations) 
			} else  { #boxcox far from 0 and 1
				result$data$obsBC <- 
						((observations^result$param["boxcox"]) - 1)/
						result$param["boxcox"]
			}
		}
	}
			


	
	
	result$model = list(reml=reml)
	if(class(trend)=="formula") {
		result$model$formula = trend
	} else {
		result$model$formula= names(trend)
	}

	
	parameterTable = data.frame(estimate=result$param)
	rownames(parameterTable) =  names(result$param)

	parameterTable$stdErr = NA
	
	stdErr = sqrt(diag(attributes(fromLogLik)$varBetaHat))
	# sometimes varBetaHat doesn't have names
	parameterTable[names(attributes(fromLogLik)$betaHat), "stdErr"] =
			stdErr

	thelims = c(0.005, 0.025, 0.05, 0.1)
	thelims = c(rbind(thelims, 1-thelims))
	
	theQ = qnorm(thelims)
	toadd = outer(parameterTable$stdErr, theQ)
	toadd = toadd + matrix(parameterTable$estimate, 
			ncol=length(thelims), nrow=dim(parameterTable)[1])
	colnames(toadd)= paste("ci", thelims, sep="")
	parameterTable = cbind(parameterTable, toadd)
	
	parameterTable[,"pval"] = pchisq(
			parameterTable$estimate^2  / parameterTable$stdErr^2,
			df=1,lower.tail=FALSE)
	
	parameterTable[,"Estimated"] = FALSE
	parameterTable[paramToEstimate,"Estimated"] = TRUE
	parameterTable[names(attributes(fromLogLik)$betaHat),"Estimated"] = TRUE
	if(estimateVariance)
		parameterTable["variance","Estimated"] = TRUE
	
	
	rownames(parameterTable)=gsub("^variance$", "sdSpatial", 
			rownames(parameterTable))	
	rownames(parameterTable)=gsub("^nugget$", "sdNugget", 
			rownames(parameterTable))	
	parameterTable[c("sdSpatial", "sdNugget"),"estimate"] = 
			sqrt(parameterTable[c("sdSpatial", "sdNugget"),"estimate"])
	
#	dimnames(parameterTable) = unlist(lapply(dimnames(parameterTable),
#			function(qq) {
#				qq=gsub("_", "\\\\textunderscore ", qq)
#				qq=gsub("\\$", "\\\\textdollar ", qq)
#				qq=gsub("<", "\\\\textless ", qq)
#				qq=gsub(">", "\\\\textgreater ", qq)
#				qq
#			}
#	))
	
	
	result$summary = as.data.frame(parameterTable)
	
	
	
	result
}


