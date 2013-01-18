krige = function(obj.model, geodata,  locations, covariates, locations.mean=locations,
		factor.info=NULL, exp.pred=FALSE,rasterMethod = c("ngb", "bilinear")) {
	
data = geodata@data	
theCovs = attributes(terms(obj.model$formula))$term.labels

varTypes = unlist(lapply(data, class))[theCovs]
thefactors = names(varTypes[varTypes=="factor"])


# if necessary, turn covariates into a raster stack
if(class(covariates) == "list") {
	
	themethod = rep(rasterMethod[1], length(covariates))
	names(themethod) = names(covariates)
	themethod[thefactors] = "ngb"
	
	locations.mean = stackRasterList(covariates,locations.mean, themethod)	
} else { # covariates must be a raster stack, use as-is
	locations.mean = covariates
}
if(!all(names(locations.mean)%in% theCovs))
	warning("some covariates in the model formula weren't supplied")

# data frame of random field prediction locations 
locationsDF = as.data.frame(locations, xy=TRUE)	

# convert data to factors, using same levels for model fit data and prediction data


# construct the fixed effects component
meanRaster = raster(locations.mean)
meanRaster[] = obj.model$beta["(Intercept)"]
for(D in theCovs){
	if(varTypes[D] == "factor") {
		if (D %in% names(factor.info)) {
			tofac = factor.info[[D]]
			if(all(c("levels","labels") %in% names(tofac)))
				tofac=data.frame(tofac[["levels"]],tofac[["labels"]])
			rownames(tofac) = tofac[,2]
			thelevels = tofac[levels(data[[D]]),1]
		}  else {
			thelevels = levels(data[[D]])	
		}
		theBetaNames = paste(D, levels(data[[D]]), sep="")
		betaHasPar = theBetaNames %in% names(obj.model$beta)
		betasHere = rep(0, length(theBetaNames))
		names(betasHere) = theBetaNames
		betasHere[betaHasPar] = obj.model$beta[theBetaNames[betaHasPar]]
		names(betasHere) = thelevels
		
		
		toAdd = setValues(locations.mean[[D]],
				betasHere[as.character(
								locations.mean[[D]]@data@values
						)]
		)
		
		meanRaster = meanRaster + toAdd
		
	} else {
		meanRaster = meanRaster + obj.model$beta[D]*locations.mean[[D]]
	}
}
names(meanRaster) = "fixed"

# do the kriging

data.col = strsplit(as.character(obj.model$formula), "~")[[2]]
geodataForKrige = as.geodata(geodata, 
		data.col=data.col, 
		covar.col=theCovs)

dummyDF = locationsDF
for(D in theCovs){
	if(varTypes[D] == "factor") {
		dummyDF[,D] = factor(levels(data[,D])[1], levels=levels(data[,D]))
	}
	else
		dummyDF[,D] = 0	
}
trend.d= trend.spatial(obj.model$trend, data)	
trend.l = trend.spatial(obj.model$trend, dummyDF)


thecontrol = krige.control(obj.model=obj.model,
		trend.d=trend.d, 
		trend.l=trend.l) 


krigeResult = krige.conv(
		geodataForKrige, locations=	locationsDF[,1:2],
		krige=thecontrol
)

rastKrige = setValues(locations[[1]], krigeResult$predict - obj.model$beta["(Intercept)"])
rastKrige = addLayer(rastKrige, 
		setValues(locations[[1]], krigeResult$krige.var)
)
names(rastKrige) = c("predict","krige.var")


if(as(meanRaster, "BasicRaster")==as(rastKrige, "BasicRaster")) {
	result = addLayer(meanRaster,
			rastKrige
	)
}	 else {
	result = addLayer(meanRaster,
			raster::resample(rastKrige, meanRaster, method="ngb")
	)
}

names(result)[names(result)=="predict"] = "random"
result = addLayer(result, 
		predict=result[["fixed"]] + result[["random"]]
)
names(result)[names(result)=="layer"] = "predict"


if(exp.pred){
	names(result)[names(result)=="predict"] = "predict.log"
	result = addLayer(result, 
			exp(result[["predict.log"]]+ 0.5*result[["krige.var"]])
	)
	names(result)[names(result)=="layer"] = "predict"
	
}

if(as(result, 'BasicRaster')!=as(rastKrige, "BasicRaster")) {
	result = list(random = rastKrige, prediction=result)			
} 		

result
}