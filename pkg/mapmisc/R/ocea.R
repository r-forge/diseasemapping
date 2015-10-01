ocea = function(x, angle=0, flip=FALSE) {

	if(!requireNamespace('geosphere', quietly=TRUE) | 
			!requireNamespace('rgdal', quietly=TRUE)) {	
		warning("install geosphere and rgdal and to use ocea")
	}	
	
	xshift = 0
	
	crsSphere = CRS(paste(
					"+proj=longlat",
					" +ellps=WGS84 +datum=WGS84", 
					" +no_defs",
					" +towgs84=0,0,0,",  
					0, # up shift
					',',
					xshift, # clockwise shift
					',',
					0, # right shift 
					',0',sep='')
	)
	
	
	if(class(x)=='Extent'){
		xExtent = projectExtent(raster(x, crs=mapmisc::crsLL), crsSphere)
	} else {
		xExtent = projectExtent(x, crsSphere)
	}
	
	midY = mean(c(ymin(xExtent),ymax(xExtent)))
	midX = mean(c(xmin(xExtent),xmax(xExtent)))
	
	myCircle = geosphere::greatCircleBearing(
			c(midX, midY),
			angle, n=5
	)
	
	myEquator = geosphere::gcLon(myCircle[2,], myCircle[3,], 0)
	
	angleIntersection=geosphere::finalBearing(
			myCircle[2,], cbind(as.vector(myEquator),0)
	)
	
	whichPos = which.max(angleIntersection)
	
	myCrs = CRS(paste("+proj=ocea",
					" +lonc=", myEquator[whichPos], 
					" +alpha=", 180-angleIntersection[whichPos], 
					" +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84", 
					" +units=m +no_defs",
					" +towgs84=0,0,0,",  
					0, # up shift
					',',
					xshift, # clockwise shift
					',',
					0, # right shift 
					',0',sep='')
	)
	
	if(flip){
		if(midY < 0) {
			myCrs = CRS(paste(as.character(myCrs), "+axis=seu"))
		} else {
			myCrs = CRS(paste(as.character(myCrs), "+axis=nwu"))
		}
	}
	
	attributes(myCrs)$crop = llCropBox(myCrs)
	
	circleLLp = SpatialPoints(
			geosphere::greatCircle(
					myCircle[2,],
					myCircle[3,], n=500, sp=FALSE
			), proj4string=crsSphere)
	
	attributes(myCrs)$circleLL = spTransform(
			circleLLp, mapmisc::crsLL
	)
	
	attributes(myCrs)$circleTrans = spTransform(
			circleLLp, myCrs)
	
	myCrs
}

moll = function(x, northShift=0, eastShift=0, twistShift=0) {
	
	
	
	northShiftS = -60*60*northShift
	eastShiftS = 60*60*eastShift
	twistShiftS = 60*60*twistShift
	
	if(any(c(northShiftS, eastShiftS, twistShiftS)!= 0)){
		resSphere	= CRS(paste(
						"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0,",
						northShiftS, ",",
						twistShiftS, ",",
						eastShiftS, ",0",
						sep=''))
	} else {
		resSphere = mapmisc::crsLL
	}
	
	if(class(x)=='Extent'){
		xExtent = projectExtent(raster(x,crs=mapmisc::crsLL))
	} else {
		xExtent = projectExtent(x, resSphere)
	}
	
	midX = mean(c(xmin(xExtent),xmax(xExtent)))
	
	res = CRS(paste("+proj=moll +lon_wrap=",
					midX, " +lon_0=",
					midX,
					" +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 ",
					"+units=m +no_defs +towgs84=0,0,0,",
					northShiftS, ",",
					twistShiftS, ",",
					eastShiftS, ",0",
					sep=''))

	
	attributes(res)$crop = llCropBox(res, 
			crsSphere=resSphere, 
			keepInner=TRUE)
	
	res
}

