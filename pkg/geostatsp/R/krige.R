

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
		forExp=as.vector(Matrix::crossprod(cholVarDataInvCovDataPred, 
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
	 	# depricated
	 	warning("no longer supported")

    
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

