ocea = function(x, angle=0, flip=FALSE, northShift=0, eastShift=0, twistShift=0) {

	if(!requireNamespace('geosphere', quietly=TRUE) | 
			!requireNamespace('rgdal', quietly=TRUE)) {	
		warning("install geosphere and rgdal and to use ocea")
	}	
	
	
	
	northShiftS = -60*60*northShift
	eastShiftS = 60*60*eastShift
	twistShiftS = 60*60*twistShift
	
	if(any(c(northShiftS, eastShiftS, twistShiftS)!= 0)){
		crsSphere	= CRS(paste(
						"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0,",
						northShiftS, ",",
						twistShiftS, ",",
						eastShiftS, ",0",
						sep=''))
	} else {
		crsSphere = mapmisc::crsLL
	}

	if(is.numeric(x)){
		if(length(x)==2) x = c(x,x+10^(-4))
		x = extent(x)
	}
	if(class(x)=='Extent'){
		xExtent = projectExtent(raster(x, crs=mapmisc::crsLL), crsSphere)
	} else {
		if(length(grep("^SpatialPoints", class(x)))){
			if(length(x)==1){
				x = raster(extend(extent(x), 10^(-4)), crs=crs(x))
			}
		}
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
					northShiftS, ",",
					twistShiftS, ",",
					eastShiftS, ",0",
					sep='') 
	)
	
	if(flip){
		if(midY < 0) {
			myCrs = CRS(paste(as.character(myCrs), "+axis=seu"))
		} else {
			myCrs = CRS(paste(as.character(myCrs), "+axis=nwu"))
		}
	}
	
	attributes(myCrs)$crop = llCropBox(
			crs=myCrs,
			crsSphere=crsSphere, keepInner=FALSE,
			N=200)
	
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


objFun = function(lonlat, target){
	
	crsRot = CRS(paste(
					"+proj=ob_tran +o_proj=longlat +o_lon_p=",
					lonlat['lon'], " +o_lat_p=",
					lonlat['lat'], 
					" +lon_0=0 +ellps=WGS84 +no_defs", sep=""))
	
	sum(as.vector(spTransform(
							target, 
							crsRot				
					)@coords*360/(2*pi))^2)
}



moll = function(x, oblique=FALSE, flip=FALSE) {
	
	
	if(is.numeric(x)){
		midX = x[1]
		if(length(x)==1) {
			midY = 0
		} else {
			midY = x[2]
		}
	} else {
		
		if(class(x)=='Extent'){
			xExtent = raster(x,crs=mapmisc::crsLL)
		} else if(length(grep("^SpatialPoints", class(x)))){
			if(length(x)==1){
				x = raster(extend(extent(x), 10^(-4)), crs=crs(x))
			}
		}
		xExtent = projectExtent(x, crsLL)
		midX = mean(c(xmin(xExtent),xmax(xExtent)))
		midY = mean(c(ymin(xExtent),ymax(xExtent)))
	}

	if(abs(midY) < 0.1 | ! oblique) {
		res = CRS(paste("+proj=moll +lon_wrap=",
					180-midX, " +lon_0=",
					midX,
					" +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 ",
					"+units=m +no_defs +towgs84=0,0,0,0,0,0,0",
					sep=''))
	} else {

		projlonlat = optim(
				c(lon=0, lat=0),
				objFun,
				target = SpatialPoints(
						cbind(midX, midY),
						proj4string=crsLL
						)
				)$par
			
		
		res = CRS(paste("+proj=ob_tran +o_proj=moll +o_lon_p=",
						projlonlat['lon'], " +o_lat_p=", projlonlat['lat'],
						" +lon_0=0 +ellps=WGS84 +datum=WGS84 ",
						"+units=m +no_defs +towgs84=0,0,0",
						sep=''))
	
	}
	
	theBox = llCropBox(
			res, 
			keepInner=FALSE, res=1)

	if(is.character(flip)) {
		res = CRS(paste(as.character(res), " +axis=", flip, sep=''))
	} else {
		if(flip){
			if(midX < 0) {
				res = CRS(paste(as.character(res), "+axis=seu"))
			} else {
				res = CRS(paste(as.character(res), "+axis=wsu"))
			}
		} 
	}
	attributes(res)$crop = theBox
	
	res
}

