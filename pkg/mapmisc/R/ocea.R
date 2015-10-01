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

moll = function(x) {
	
	
	if(class(x)=='Extent'){
		xExtent = x
	} else {
		xExtent = projectExtent(x, mapmisc::crsLL)
	}
	
	midX = mean(c(xmin(xExtent),xmax(xExtent)))
	
	res = CRS(paste("+proj=moll +lon_wrap=",
					midX, " +lon_0=",
					midX,
					" +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs84=0,0,0",
					sep=''))
	
	attributes(res)$crop = llCropBox(res)
	
	res
}

