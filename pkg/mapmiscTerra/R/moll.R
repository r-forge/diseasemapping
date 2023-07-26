

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
	
	# x = vect(cbind(-100,45), crs=crsLL)

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
			if(terra::is.points(x) & length(x)==1){
				x = rast(terra::extend(terra::ext(x), 10^(-4)), crs=crs(x))
			}
		}
		xExtent = project(terra::ext(x), terra::crs(x), crsLL)
		midX = mean(c(terra::xmin(xExtent),terra::xmax(xExtent)))
		midY = mean(c(terra::ymin(xExtent),terra::ymax(xExtent)))
	}
	xVec = vect(cbind(midX, midY), crs = crsLL)

	if(is.null(angle)){
		angle = 0
		result = paste(
			"+proj=moll +lon_wrap=",
			midX, " +lon_0=",
			midX,
			" +x_0=0 +y_0=0 ",
			"+datum=WGS84 +ellps=WGS84 ",
			"+units=m +no_defs",
			sep='')
	} else {

		testPoints = vect(c(xVec, 
			vect(geosphere::destPoint(terra::crds(xVec), b=angle, d=500*1000), 
				crs=crs(xVec))))

		testFunCrs = function(param, testPoints) {
			crsHere = crsFromNumeric(param)
			xxx = try(terra::crds(project(testPoints, crsHere)), silent=TRUE)
			if(all(class(xxx) == 'try-error')) xxx = matrix(NA, 2,2)
				xxx = xxx/1e6
	# x coords the same, x point close to origin, dest point above x
			sqrt(diff(xxx[,1])^2 + sum(xxx[1, ]^2) + abs(xxx[2,1]-xxx[1,1])^2)
		}
		crsFromNumeric = function(param) {

			paste(
				'+proj=ob_tran +o_proj=moll ',
				'+o_lon_p=',  param['lon_p'],	
				'+o_lat_p=', param['lat_p'],  
				'+lon_0=', sum(param[c('lon_0minuslon_p', 'lon_p')])
			)

		}


		fromOptim = optim(
			c(lon_p=midX,lat_p=midY, lon_0minuslon_p=0),
			testFunCrs,
			lower = c(-180, -90, -100), upper = c(180, 90, 100),
			method = 'L-BFGS-B',
			testPoints = testPoints
		)

		result = crsFromNumeric(fromOptim$par)				
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

	phiSeq = seq(0, 2*pi, len=1001)
	theBox=list(ellipse=vect(cbind(1.8e7*cos(phiSeq), 9e6*sin(phiSeq)), type='polygons', crs=result))
	if(angle ==0) {
		angleOpposite = Arg(exp(1i*2*pi*(midX + 180)/360))*(360/(2*pi))
		densify.interval = 10*1000
		theLine = terra::densify(vect(rbind(cbind(angleOpposite,-89.9), cbind(angleOpposite, 89.9)), type='lines', crs=crsLL), densify.interval)
		theBox$crop = terra::buffer(theLine, width=50*1000)
		theBox$regionLL = terra::erase(terra::as.polygons(terra::ext(-100,180,-90,90), crsLL), theBox$crop)
	} else {
		theBox = llCropBox(result, ellipse = theBox$ellipse, 
			buffer.width=75*1000, densify.interval = 20*1000, 
  			crop.distance = Inf, crop.poles = FALSE, crop.leftright=FALSE,
		  remove.holes=TRUE)
	}

	# to do: ellipse wrong when flip is used

	attributes(result)[names(theBox)] = theBox

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
		attributes(myCrs)$regionLL = cropBox$regionLL


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
