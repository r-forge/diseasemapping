

# y is gaussian, w is positive
# y = (w^lambda - 1)/lambda
# lambda * y + 1 = w^lambda = exp(lambda * log(w))

meanBoxCox=function(pred, sd, boxcox, Nbc = 100) {
  
  # if boxcox is negative, don't let values get closer
  # than this to zero
  epsX = exp(12*boxcox)
  
  epsBC = 0.001 # check box-cox within this distance
  # of zero or one
  if(abs(boxcox)<epsBC){
    return(
      exp(pred  + sd^2/2)
    )
  } else if (abs(boxcox-1)<epsBC){
    return(pred)
  }
  
  # normal probabilities for numerical integration
  SXboxcox = seq(-7,7,len=Nbc)
  SXboxcox = unique(signif(
      c(
        SXboxcox,
        SXboxcox/SXboxcox[1]
      ),6)
  )
  SXboxcox = sort(SXboxcox)
  PXboxcox = pnorm(SXboxcox)#, log=TRUE)
  # probability of normal being in a bin
  DXboxcox = diff(PXboxcox)/2
  NDX = length(DXboxcox)
  DXboxcox = c(
    DXboxcox[1],
    DXboxcox[-1] + DXboxcox[-NDX],
    DXboxcox[NDX]
  )
  IXboxcox = log(DXboxcox)  
  
  x = boxcox * (outer(sd, SXboxcox) +  pred)+1
  
  # negatives to zero
  xneg= which( as.vector(x < 0))
  
  if(boxcox<0){
    # get rid of values very close to zero
    xneg = c(xneg,
      which(abs(as.vector(x)) < epsX)
    )
  }
  
  x[xneg] = NA
  
  logx = log(x)/boxcox 
  
  IXmat = matrix(IXboxcox, nrow(x), ncol(x), byrow=TRUE)
  
  result =  rowSums(exp(logx + IXmat),na.rm=TRUE)
  allNA = rowSums(!is.na(x))==0
  
  if(length(xneg)) {
    IXmat[xneg] = NA
    result = cbind(
      predict=result,
      probComplex.boxcox = 
        1-rowSums(exp(IXmat),na.rm=TRUE)
    )
    result[allNA,] = NA        
  } else {
    result[allNA] = NA
  }
  
  
  result
}


krigeOneRowPar = function(
	Drow, yFromRowDrow, 
	locations,
	param,coordinates,Ny,
	cholVarDataInv,
	cholVarDatInvData,
	xminl,xresl,ncoll,
	lengthc){

		  # covariance of cells in row Drow with data points
	resC =  .C(C_maternArasterBpoints, 
		as.double(xminl), 
		as.double(xresl), 
		as.integer(ncoll), 
		as.double(yFromRowDrow), 
		as.double(0), as.integer(1),
		as.double(crds(coordinates)[,1]), 
		as.double(crds(coordinates)[,2]), 
		N=as.integer(Ny), 
		result=as.double(matrix(0, ncoll, 
			lengthc)),
		as.double(param["range"]),
		as.double(param["shape"]),
		as.double(param["variance"]),
		as.double(param["anisoRatio"]),
		as.double(param["anisoAngleRadians"])
	) 
	covDataPred = matrix(resC$result, nrow=ncoll, ncol=Ny)


	cholVarDataInvCovDataPred = tcrossprod(cholVarDataInv, covDataPred)

	x= cbind( # the conditional expectation
		forExp=as.vector(crossprod(cholVarDataInvCovDataPred, 
			cholVarDatInvData)),
				  # part of the conditional variance
		forVar=apply(cholVarDataInvCovDataPred^2, 2, sum)
	) 
	x

}


krigeLgm = function(
		formula, data, 
		grid,
		covariates=NULL, 
		param,   
		expPred=FALSE,
		nuggetInPrediction=TRUE, 
		mc.cores=getOption("mc.cores", 1L)) {
  
	 # this function really needs some tidying!
	 
	 trend = formula
	 locations = grid
	 coordinates=data
	 theVars = NULL
	 
	 haveBoxCox = any(names(param)=="boxcox")
	 NsimBoxCox=50
  if(haveBoxCox) {
    haveBoxCox = abs(param["boxcox"]-1) > 0.001
    if(param['boxcox']<0) NsimBoxCox=100
    if(param['boxcox']< -0.2) NsimBoxCox=200
    if(param['boxcox']< -0.5) NsimBoxCox=400
  }
	 
	 haveNugget = any(names(param)=="nugget")
	 if(haveNugget) { 
		  haveNugget = param["nugget"] > 0
	 } 
	 if(!haveNugget) {
		  nuggetInPrediction=FALSE
	 }
	 
	 if(is.numeric(locations)){
		  locations = squareRaster(data, locations)
	 }
	 if(nrow(locations) * ncol(locations) > 10^7) warning("there are lots of cells in the prediction raster,\n this might take a very long time")
	 
	 
	 observations = meanRaster = NULL
	 
	 noCovariates = length(names(covariates)) == 0
	 if(is.data.frame(covariates)) {
	 	if(!nrow(covariates)) noCovariates = TRUE
	 }
  
	 if(noCovariates) {
		  # no coariates, mean is intercept
		  if(any(names(param)=='(Intercept)')) {
			   meanForRaster = param['(Intercept)']		
		  } else {
			   meanForRaster = 0
		  }
		  meanFixedEffects = 
				  rep(meanForRaster, ncell(locations))
		  meanRaster = terra::rast(locations)
		  terra::values(meanRaster) = meanFixedEffects	
	 }

	 if( is.data.frame(covariates) & any(class(formula)=="formula") & !noCovariates)  {
    
		    # put zeros for covariates not included in the data frame
		    notInCov = setdiff(all.vars(formula), names(covariates))
		    for(D in notInCov)
			     covariates[[D]] = 0
      
		    modelMatrixForRaster = model.matrix(formula, covariates)
		    if(nrow(modelMatrixForRaster) < nrow(covariates)) {
		    	# some cells missing covariates
		    	toAdd = setdiff(rownames(covariates), rownames(modelMatrixForRaster))
		    	toAdd = matrix(NA, length(toAdd), ncol(modelMatrixForRaster),
		    			dimnames = list(toAdd, colnames(modelMatrixForRaster)))
		    	modelMatrixForRaster = rbind(
		    		modelMatrixForRaster, toAdd
		    	)[rownames(covariates), ]
		    }
      
		    theParams = intersect(colnames(modelMatrixForRaster), names(param))
		    
		    meanForRaster = drop(
				    tcrossprod( param[theParams], modelMatrixForRaster[,theParams] )
		    )
		    meanFixedEffects = rep(NA, ncell(locations))
		    if('space' %in% names(covariates)  & 'space' %in% names(locations)) {
		    	meanFixedEffects[match(covariates$space, terra::values(locations)[,'space'])] = meanForRaster
		    } else {
			    meanFixedEffects[as.integer(names(meanForRaster))] = meanForRaster
			  }
		    meanRaster = rast(locations)
		    terra::values(meanRaster) = meanFixedEffects
      
      
	 } # end covariates is DF
	 
	 if(any(class(data)=="SpatVector")&
	 		any(class(formula)=="formula")) {
	   
		  if(all(setdiff(names(covariates), 'space') %in% names(data))) {
      
			   modelMatrixForData = model.matrix(formula, values(data))
      
			   theParams = intersect(colnames(modelMatrixForData), names(param))
      
			   meanForData =	 as.vector(tcrossprod(
						    param[theParams],
						    modelMatrixForData[,theParams])
			   )
			   names(meanForData) = rownames(modelMatrixForData)
			   
			   haveData = match(names(meanForData), 
				    rownames(values(data)))
      
      data = data[haveData,]
			   coordinates=data
      
			   
		    observations = drop(values(data)[,
						    all.vars(formula)[1] ] )
		    
		    if(haveBoxCox) {
			     if(abs(param["boxcox"]) < 0.001) {
				      observations = log(observations)
				      expPred = TRUE
				      haveBoxCox = FALSE
			     } else {
				      observations = ((observations^param["boxcox"]) - 1)/
						      param["boxcox"]
			     }
			     
		    } # have boxcox
		    observations = observations - meanForData		
		  } # end all covariates in data
	   
	 } # end data is spdf	
	 
	 
	 if(!length(observations) | is.null(meanRaster)) {


	  # old code, not called from lgm
		  # the above didn't create observations and meanRaster
		  # use the old code, probably not being called from lgm
 	  
	   # find factors, so we reproject rasters using
	   # the correct method.
	   # search for factors in the data supplied
    
    
	   
	   # look for factors in the model formula
	   if(any(class(trend)=="formula")){
      
      trendFormula = update.formula(trend, junk ~ . )
      
      
		    covariatesForData = values(data)
 		   
		    if(is.vector(data)) {
			     observations = data
		    } else {			
			     observations = all.vars(trend)[1]
			     observations = covariatesForData[,observations]
		    }
      
		    theVars = all.vars(trendFormula)[-1]
      
		    if(length(theVars)) {
			     factorsInData = unlist(lapply(
							     covariatesForData[,theVars,drop=FALSE],
							     is.factor))
			     factorsInData = names(factorsInData)[factorsInData]
		    } else {
			     factorsInData=NULL
		    } 
		    
		    
		    allterms = rownames(attributes(terms(trend))$factors)
		    
		    factorsInFormula = grep("^factor", allterms, value=TRUE)
		    factorsInFormula = gsub("^factor\\(", "", factorsInFormula)
		    factorsInFormula = gsub("\\)$", "", factorsInFormula)
		    
		    factorsInTrend=NULL
		    
		    allterms = gsub("^[[:alnum:]]+\\(", "", allterms)
		    allterms = gsub("\\)$", "", allterms)
		    
		    if(!all(allterms %in% names(data)))
			     warning("some covariates don't appear in data")
	   } else { # trend not formula
		    # trend is a data frame of covariates
		    # look for factors in it
		    covariatesForData = as.data.frame(trend)
		    
		    observations = as.data.frame(data)[,1]
		    
		    factorsInTrend = unlist(lapply(
						    covariatesForData, is.factor
				    ))
		    factorsInTrend = names(factorsInTrend)[factorsInTrend]
		    factorsInFormula = factorsInData = NULL
		    
		    # guess at the formula
		    trendFormula = as.formula(paste(
          "junk ~ ",
						    paste(c('1', names(covariatesForData)), collapse="+")
				    )
      )
	   } # end trend not a formula
	   
   
	   # we know which variables factors
	   theFactors = unique(c(factorsInFormula, factorsInData, factorsInTrend))
	   theFactors = theFactors[theFactors %in% names(covariates) ]
    
 	  if(length(grep("SpatRaster|^list", class(covariates)))) { 
 	    # if there's only variable in the model assign it's name to covariates
	     covariateNames = all.vars(
        update.formula(trendFormula, junk~ . )
      )[-1]
	     if(length(covariateNames)==1){
		      # so far only one variable
		      names(covariates)= covariateNames
	     } 
	     # loop through factors
	     # and make sure integer values in rasters get converted
	     # to things with parameter values!
	     for(D in theFactors) {
		      # is this variable in param with  a factor around it? 
		      # for instance factor(x)1 and factor(x)2 ?
		      paramWithFactor = grep(
				      paste("factor\\(", D, "\\)[[:digit:]]+$", sep=""),
				      names(param), value=TRUE)
		      paramStartWithD = grep(
				      paste("^", D, ".+$", sep=""),
				      names(param), value=TRUE)
		      paramFactorCharacter = grep(
				      paste("factor\\(", D, "\\).+$", sep=""),
				      names(param), value=TRUE)
		      if(length(paramWithFactor)) {
			       # formula will convert to factor, don't 
			       # create factor beforehand
			       theLevels = gsub(
					       paste("^factor\\(",D,"\\)",sep=""),
					       "",paramWithFactor)
			       theLevels = as.integer(theLevels)
			       allValues = unique(covariates[[D]])[,1]
			       dontHave = allValues[!allValues %in% theLevels]
			       # make values with no data all equal to the lowest value
		        # so it's the baseline when turning into a factor.
			       forRecla = cbind(dontHave, min(allValues)-1)
			       
			       covariates[[D]] = 
					       classify(covariates[[D]], forRecla)
	         
		      } else if( length(paramStartWithD) ) {
			       # not a bunch of digits, 
			       # stuff like xTrees and xGrassland for covariate x and levels Trees and Grassland
			       # see if these line up with 
			       theLevels = gsub(paste("^", D, sep=""),"",paramStartWithD)

			       levelsTable = levels(covariates[[D]])[[1]]

			       inId = theLevels %in% as.character(levelsTable[,1])
			       inLabel = theLevels %in% levelsTable[,2]
          
			       if(mean(inId) > mean(inLabel)){
				        levelsTable$levelsInParams =  
						        as.character(levelsTable[,1])
				        labelCol = ncol(levelsTable)
				        levelsInTable = levelsTable[,1] %in% 
						        theLevels
			       } else {
				        levelsInTable = levelsTable[,2]%in% theLevels
				        labelCol=2
			       }
			       
			       if(mean(theLevels %in% levelsTable[,labelCol]) < 0.4)
				        warning("many levels appear missing in covariate", D)
			       valuesInParams = levelsTable[levelsInTable,1]
          
			       allValues = unlist(unique(covariates[[D]]))
			       allValues = levelsTable[levelsTable[,2] %in% allValues, 1]
			       dontHave = allValues[!allValues %in% valuesInParams]
			       forRecla = cbind(dontHave, min(as.vector(minmax(covariates[[D]])))-1)
			       covariates[[D]] = classify(covariates[[D]], forRecla)
			       
			       
			       levelsTable = 
					       levelsTable[c(1, 1:nrow(levelsTable)),c(1,labelCol)]
			       levelsTable[1,1]= min(allValues)-1
			       levelsTable[1,2] = ''
			       colnames(levelsTable)[2] = "levels"
			       levels(covariates[[D]]) = levelsTable
			       
		      } else if (length(paramFactorCharacter)) {
			       # stuff like factor(x)Trees and factor(x)Grassland for covariate x and levels Trees and Grassland
			       theLevels = gsub(paste("^factor\\(", D,"\\)", sep=""),"",
					       paramFactorCharacter)
			       levelsTable = levels(covariates[[D]])[[1]]
			       levelsInTable = levelsTable[,2]%in% theLevels
			       if(mean(theLevels %in% levelsTable[,2]) < 0.4)
				        warning("many levels appear missing in covariate", D)
			       valuesInParams = as.numeric(levelsTable[levelsInTable,1])
          
#			       allValues = unlist(unique(covariates[[D]]))
#			       dontHave = allValues[!allValues %in% valuesInParams]
			       forRecla = cbind(levelsTable[!levelsInTable,'ID'], min(levelsTable[,'ID'])-1)
			       covariates[[D]] = classify(covariates[[D]], forRecla)
			       
			       
			       levelsTable = 
					       levelsTable[c(1, 1:nrow(levelsTable)),]
			       levelsTable[1,1]= min(levelsTable[,'ID'])-1
			       levelsTable[1,2] = "0"
			       colnames(levelsTable)[2]="levels"
			       levels(covariates[[D]])[[1]] =  levelsTable			
			       
			       
		      } else {
			       warning("don't know what to do with covariate", D, 
					       "\n can't assign parameters to levels of this factor")			
		      }
		      
		      
	     } # end loop through factors
	     
	     if(length(grep("SpatRaster|^list", class(covariates))) & length(theVars)) {
		      # method for resampling covariate rasters
        
		      method = resampleMethods(formula, covariates)
		      
		      covariates = stackRasterList(covariates, template=locations, method=method)
        
		      theVars = do.call('intersect',
          dimnames(attributes(terms(trendFormula))$factors))
        
		      if(nlyr(covariates)==1 & length(theVars)==1) {
			       names(covariates) = theVars
		      }
		      
		      # construct the fixed effects component
		      covariatesDF = cbind(values(covariates, dataframe=TRUE), 
		      	crds(covariates, df=TRUE, na.rm=FALSE))
		      # get rid of trailing _ created by as.data.frame
		      names(covariatesDF) = gsub("_levels$", "", names(covariatesDF))
	     } else {
		      covariatesDF = as.data.frame(matrix(NA, ncol=0, nrow=ncell(locations)))
	     }
	   } else {# end covariates is raster, assume it's a data frame
		    covariatesDF=covariates
	   } 
    
	   # get rid of response variable in trend formula
    meanRaster = rast(locations)
    names(meanRaster) = "fixed"
    
	   
	   if(length(all.vars(trendFormula)) ){ # if have covariates
	     missingVars = all.vars(trendFormula)[-1] %in% names(covariatesDF)
	     missingVars = all.vars(trendFormula)[-1][!missingVars]
	     
	     # check if all variables are in covariates
	     if(length(missingVars)) {
		      cat("cant find covariates ",
				      paste(missingVars, collapse=","),
				      " for prediction, imputing zeros\n")		
		      
		      covariatesDF[,missingVars]=0	
	     }
      
	     modelMatrixForRaster = model.matrix(trendFormula, cbind(covariatesDF,junk=0))
      
	     if(!all(colnames(modelMatrixForRaster)%in% names(param))){
		      warning("cant find coefficients",
				      paste(setdiff(colnames(modelMatrixForRaster), names(param)), collapse=", "),
				      "in param\n")
	     }
	     
	     
	     meanFixedEffects = 
			     modelMatrixForRaster %*% param[colnames(modelMatrixForRaster)]
	     
	     anyNA = apply(covariatesDF, 1, function(qq) any(is.na(qq)))
	     if(any(anyNA)) {
		      oldmm = rep(NA, ncell(meanRaster))
		      oldmm[!anyNA] = meanFixedEffects
		      terra::values(meanRaster) = oldmm
	     } else {
		      terra::values(meanRaster) = meanFixedEffects
	     }
	     
      
	     
	     modelMatrixForData = model.matrix(trendFormula, 
        cbind(covariatesForData,junk=0))
	     haveData = match(rownames(modelMatrixForData), 
			     rownames(covariatesForData))
	     observations = observations[haveData]
	     coordinates = coordinates[haveData,]
	     
	     meanForData = 
			     modelMatrixForData %*% param[colnames(modelMatrixForData)]
      
	   } else { #no covariates	
      
		    if(any(names(param)=='(Intercept)')) {
			     terra::values(meanRaster) = param['(Intercept)'] 
		    } else {
			     terra::values(meanRaster) = 0
		    }
		    meanForData = rep(values(meanRaster)[1], length(observations))
	   }
	   
    
	   
# subtract mean from data
    
	   theNAdata =  
			   is.na(observations)
	   
	   if(all(theNAdata)) {
		    warning(
				    'it appears there are no observations without at least one covariate missing')
	   }
		  
    
	   
	   if(any(theNAdata)) {
		    noNAdata = !theNAdata
		    if(length(grep("^SpatVector", class(coordinates)))) {
			     coordinates = coordinates[noNAdata,]	
		    } else if(any(class(coordinates)=="dist")){
			     coordinates = as.matrix(coordinates)
			     coordinates = coordinates[noNAdata,noNAdata]
			     coordinates = as.dist(coordinates)
		    } else {
			     warning("missing vlaues in data but unclear how to remove them from coordinates")
		    }
		    observations = observations[noNAdata]
	   }
	   
	   
    
	   
	   if(haveBoxCox) {
		    if(abs(param["boxcox"]) < 0.001) {
			     observations = log(observations)
			     expPred = TRUE
			     haveBoxCox = FALSE
		    } else {
			     observations = ((observations^param["boxcox"]) - 1)/
					     param["boxcox"]
		    }
		    
	   } # end have box cox
	   
	   observations = observations - meanForData
    

    
	 } # end old code not called from LGM

	 	 
  cholVarDataInv = geostatsp::matern(
  	coordinates, 
  	param=param, type='inverseCholesky')

    
  cholVarDatInvData = cholVarDataInv %*% observations

	 Ny = length(observations)
	 param = fillParam(param)


	 datForK = list(
			 locations=locations,param=param,
			 coordinates=coordinates,Ny=Ny,
			 cholVarDataInv=cholVarDataInv,
			 cholVarDatInvData = cholVarDatInvData,
			 xminl=xmin(locations),
			 xresl = xres(locations),
			 ncoll=ncol(locations),
			 lengthc=length(coordinates)
		
		)
	 Srow = 1:nrow(locations)
	 
	 if(mc.cores ==1 ) {
	   sums=mapply(krigeOneRowPar, 
	   	Drow = Srow, 
				  yFromRowDrow = yFromRow(locations,Srow),
				  MoreArgs=datForK,
				  SIMPLIFY=FALSE)
    
	 } else {
		  sums=parallel::mcmapply(krigeOneRowPar, Srow, 
				  yFromRow(locations,Srow),
				  MoreArgs=datForK,SIMPLIFY=FALSE,mc.cores=mc.cores)
		  
	 }

   
 	sums <- simplify2array(sums)
	 # row sums of cholVarDataInvCovDataPred
	 forExpected = sums[,'forExp',]
	 # row sums of squares
	 forVar = sums[,'forVar',]
	 
	 randomRaster = rast(meanRaster)
	 names(randomRaster) = "random"
 	terra::values(randomRaster) = t(forExpected)
	 

	 predRaster = meanRaster + randomRaster
	 names(predRaster) = "predict"
	 
	 if(any(forVar > param["variance"])){
		  themax = max(forVar - param["variance"],na.rm=TRUE)
		  if(themax > 1e-6)
			   warning("converted variances of ", themax, " to zero")	
		forVar = pmin(forVar, param["variance"])	
	 }
	 
	 krigeSd = rast(meanRaster)
	 names(krigeSd) = "krigeSd"
  
  if(nuggetInPrediction) {
		  terra::values(krigeSd) = sqrt(sum(param[c("nugget","variance")]) - 
						  as.vector(forVar))
	 } else {
		  terra::values(krigeSd) = sqrt(param["variance"] - as.vector(forVar))
	 }
	 
  names(meanRaster) = "fixed"

  
	 result = c(meanRaster, randomRaster, predRaster,
			 krigeSd)
  
	 # box-cox
	 if(haveBoxCox){ 
		  names(result)[names(result)=="predict"] = "predict.boxcox"
		  
    bcpred = meanBoxCox(
      pred=values(result[['predict.boxcox']], mat=FALSE, dataframe=FALSE), 
      sd= values(result[['krigeSd']], mat=FALSE, dataframe=FALSE),
      boxcox=param['boxcox']
    )
    
    newraster=rast(result[["predict.boxcox"]])
    names(newraster) = "predict"
    if(is.matrix(bcpred)){
      terra::values(newraster) = bcpred[,'predict']
      add(result) = newraster
#      names(newraster) = 'probComplex.boxcox'
#      values(newraster) = bcpred[,'probComplex.boxcox']
#      result = addLayer(result, 
#          newraster)
    } else {
      terra::values(newraster) = bcpred
      add(result) = newraster
    }
    
    
		  
	 } # end have box cox
	 
	 
	 if(expPred){
		  
		  names(result)[names(result)=="predict"] = "predict.log"
		  newLayer = exp(result[["predict.log"]]+ 0.5*result[["krigeSd"]]^2 )
		  names(newLayer) = "predict"
      add(result) = newLayer

		  
	 } # end expPred
		
	 result
}

