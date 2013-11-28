

modelRandomFields = function(param){
	

param = fillParam(param)

param["scaleRandomFields"] = param["range"]/2 

if(abs(param["anisoRatio"]) >  10^(-4)){
	# geometric anisotropy
	scaleTwoDirections = c(1/param["scaleRandomFields"],
			1/(param["anisoRatio"]*param["scaleRandomFields"]))
	angle = param["anisoAngleRadians"]
	anisoMat = diag(scaleTwoDirections) %*% 
			matrix(c(cos(angle), sin(angle), 
							-sin(angle), cos(angle)),2)
	model=list("$", var=param["variance"],   
			A=anisoMat,
			list("matern", nu=param["shape"]))	
} else {
	model=list("$", var=param["variance"],   
			s=param["scaleRandomFields"],
			list("matern", nu=param["shape"]))	
}	

model 
}

fillParam = function(param) {

	if(is.numeric(param))
		param = matrix(param, ncol=length(param), nrow=1,
				dimnames=list("1", names(param)))

	if(is.list(param)) {

		parlengths = unlist(lapply(param, length))
		parlengths = unique(parlengths)
		Nsamples = parlengths[parlengths != 1]
		if(length(Nsamples)!= 1) {
			warning("some parameters have more samples than others")
		}
	
		param = lapply(param, as.vector)
		param = do.call(cbind, param)
	}

	
	colnames(param) = gsub("^var$", "variance", colnames(param))
	
	sdname = grep("^sd",colnames(param), ignore.case=TRUE,value=TRUE)
	if(!any(colnames(param)=="variance") & 
			length(sdname)) {
		param = cbind(param, variance = param[,sdname]^2)
	}
# if still no variance set it to 1.
	if(!any(colnames(param)=="variance")) {
		param=cbind(param, variance = 1)
	}

	if(!any(colnames(param)=="nugget")) {
		if(any(colnames(param)=="nuggetSD")) {
			param = cbind(param, nugget = param[,"nuggetSD"]^2)
		}	else {
			param = cbind(param, nugget=0)
		}
	}
	# shape
	if(!any(colnames(param)=="shape"))
		warning("shape parameter not supplied")
	
# fill in anisotropy parameters	
	if(!any(colnames(param)=="anisoRatio"))
		param = cbind(param, anisoRatio = 1)
	if(!any(colnames(param)=="anisoAngleRadians")){
		if(any(names(param)=="anisoAngleDegrees")) {
			param = cbind(param, 
					anisoAngleRadians=
					param[,"anisoAngleDegrees"]*2*pi/360)
		} else {
			param = cbind(param,anisoAngleRadians = 0)
		}
	}
	drop(param)
}
