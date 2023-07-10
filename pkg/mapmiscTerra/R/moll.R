

objFun = function(param, start, target){
	
	crsRot = crs(paste(
					"+proj=ob_tran +o_proj=robin +o_lon_p=",
					param['lon'], " +o_lat_p=",
					param['lat'], 
					" +ellps=sphere +no_defs", sep=""))
	
	sum(
			(as.vector(terra::crds(project(
							start, 
							crsRot				
					)))*360/(2*pi)
			 - target
			)^2, na.rm=TRUE
	)
}


objFunAngle = function(param, start, target) {

	# start is a SpatialPoints object of length 2
	# target is x, y
	# point 1 should be due north of point 2 after transform
	
	crsRot = crs(paste("+proj=ob_tran +o_proj=moll +o_lon_p=",
					param['lon'], " +o_lat_p=", param['lat'],
					" +lon_0=", param['wrap'], 
					" +lon_wrap=", param['wrap'],
					" +ellps=WGS84 ",
					"+units=m +no_defs +towgs84=0,0,0",
					sep='')
	)
	
	startT = project(
					start, 
				crsRot				
	)
	
	resDist = sum(
			(as.vector(terra::crds(startT)[1,] - target[c('x','y')])/100000
						)^2
	)
	
	resAngle = terra::crds(startT)[2,] - terra::crds(startT)[1,]
	resAngle = 90 - Arg(resAngle[1] + 1i * resAngle[2])*360/(2*pi) 
	
	resDist + resAngle^2
	
}

moll = function(x=0, angle=NULL, flip=FALSE) {
	
	
	if(is.numeric(x)){
		midX = x[1]
		if(length(x)==1) {
			midY = 0
		} else {
			midY = x[2]
		}
	} else {
		
		if(any(class(x)=='SpatExtent')){
			xExtent = rast(x,crs=crsLL)
		} else if(length(grep("^SpatVector", class(x)))){
			if(length(x)==1){
				x = rast(terra::extend(terra::ext(x), 10^(-4)), crs=crs(x))
			}
		}
		xExtent = project(terra::ext(x), terra::crds(x), crsLL)
		midX = mean(c(terra::xmin(xExtent),terra::xmax(xExtent)))
		midY = mean(c(terra::ymin(xExtent),terra::ymax(xExtent)))
	}

	if(is.null(angle)){
		result = paste(
						"+proj=moll +lon_wrap=",
						midX, " +lon_0=",
						midX,
						" +x_0=0 +y_0=0 ",
						"+datum=WGS84 +ellps=WGS84 ",
						"+units=m +no_defs",
						sep='')
	} else {

	newPole = geosphere::destPoint(
			c(midX, midY), b=angle, 
			d=pi/2, a=1, f=0
	)
	
#	stuff=geosphere::greatCircle(c(midX, midY), newPole, n=100)
	
#	geosphere::onGreatCircle(c(midX, midY), c(0,90), newPole)
	
	param = optim(
				c(lon=0,lat=0),
				objFun,
				start = vect(
						newPole,
						crs=crsLL
						),
				target=c(NA,10^9)
				)
				
	param = param$par			
			
#	spTransform(
#			SpatialPoints(newPole, proj4string=crsLL),
#			CRS(paste(
#							"+proj=ob_tran +o_proj=longlat +o_lon_p=",
#							param['lon'], " +o_lat_p=",
#							param['lat'], 
#							" +lon_0=0 +ellps=WGS84 +no_defs", sep=""))
#	)@coords*360/(2*pi)
	
	
	newOrigin = terra::crds(project(
			vect(cbind(midX, midY), crs=crsLL),
			paste(
							"+proj=ob_tran +o_proj=longlat +o_lon_p=",
							param['lon'], " +o_lat_p=",
							param['lat'], 
							" +lon_0=0 +ellps=WGS84 +no_defs", sep="")
	))*360/(2*pi)

	# optimize again, with angle

	paramAgain = c(param, wrap=as.numeric(newOrigin[1]))

	start = cbind(midX, midY)
	
	start = vect(rbind(
					start, 
					geosphere::destPoint(start, angle, d=1000)
					), crs=crsLL)
			
	newParam = optim(paramAgain, objFunAngle, start=start, 
			target=c(x=0, y=0, angle=as.numeric(angle))
	)$par		
		
	result = paste("+proj=ob_tran +o_proj=moll +o_lon_p=",
						newParam['lon'], " +o_lat_p=", newParam['lat'],
						" +lon_0=", newParam['wrap'], 
						" +lon_wrap=", newParam['wrap'],
						" +ellps=WGS84 ",
						"+units=m +no_defs +towgs84=0,0,0",
						sep='')
	
	}
	

	if(is.character(flip)) {
	  result = paste(as.character(result), " +axis=", flip, sep='')
	} else {
		if(flip){
			if(midX < 0) {
			  result = paste(as.character(result), "+axis=seu")
			} else {
			  result = paste(as.character(result), "+axis=wsu")
			}
		} 
	}
  result  = crs(result)

  theBox = llCropBox(crs=result)
  
#	 attributes(result)$regionLL = theBox$poly
		attributes(result)$ellipse = theBox$ellipse
	 
	 attributes(result)$crop = theBox$crop
	 
#	 attributes(result)$ellipseSafeLL = theBox$polyTrans
	 
	 
	result
}

ocea = function(x, angle=0, flip=FALSE) {
	
	
	northShift=0; eastShift=0; twistShift=0
	northShiftS = -60*60*northShift
	eastShiftS = 60*60*eastShift
	twistShiftS = 60*60*twistShift
	
	if(any(c(northShiftS, eastShiftS, twistShiftS)!= 0)){
		crsSphere	= crs(paste(
						"+proj=longlat +no_defs +ellps=WGS84 +towgs84=0,0,0,",
						northShiftS, ",",
						twistShiftS, ",",
						eastShiftS, ",0",
						sep=''))
	} else {
		crsSphere = crsLL
	}
	
	if(is.numeric(x)){
		if(length(x)==1) x = c(x,0)
		if(length(x)==2) x = c(x,x+10^(-2))
		x = terra::ext(x)
	}
	if(any(class(x)=='SpatExtent')){
		xExtent = project(x, crsLL, crsSphere)
	} else {
		if(length(grep("^SpatVector", class(x)))){
			if(length(x)==1){
				x = rast(terra::extend(terra::ext(x), 10^(-4)), crs=crs(x))
			}
		}
		xExtent = project(x, crs(x), crsSphere)
	}
	
	midY = mean(c(terra::ymin(xExtent),terra::ymax(xExtent)))
	midX = mean(c(terra::xmin(xExtent),terra::xmax(xExtent)))
	
	myCircle = geosphere::greatCircleBearing(
			c(midX, midY),
			angle, n=5
	)
	if(all(is.na(myCircle[,2]))) myCircle[,2] = midY
	
	myEquator = geosphere::gcLon(myCircle[2,], myCircle[3,], 0)
	
	angleIntersection=geosphere::finalBearing(
			myCircle[2,], cbind(as.vector(myEquator),0)
	)
	
	whichPos = which.max(angleIntersection)
	
	myCrs = paste("+proj=ocea",
					" +lonc=", myEquator[whichPos], 
					" +alpha=", 180-angleIntersection[whichPos], 
					" +x_0=0 +y_0=0 +ellps=WGS84", 
					" +units=m +no_defs",
					" +towgs84=0,0,0,",  
					northShiftS, ",",
					twistShiftS, ",",
					eastShiftS, ",0",
					sep='') 

	if(identical(flip, TRUE)){
		if(midY >= 0) {
			myCrs = paste(as.character(myCrs), "+axis=esu")
		} else {
			myCrs = paste(as.character(myCrs), "+axis=nwu")
		}
	}
	if(is.character(flip))
		myCrs = paste0(as.character(myCrs), " +axis=", flip)

	myCrs = crs(myCrs)

	cropBox = llCropBox(crs=myCrs, crop.distance = Inf)
	
	attributes(myCrs)$crop = cropBox$crop
	attributes(myCrs)$ellipse = cropBox$ellipse

	
	circleLLp = vect(
			geosphere::greatCircle(
					myCircle[2,],
					myCircle[3,], n=500, sp=FALSE
			), crs=crsSphere)
	
	attributes(myCrs)$circleLL = project(
			circleLLp, crsLL
	)
	
	attributes(myCrs)$circleTrans = project(
			circleLLp, myCrs)
	
	myCrs
}
