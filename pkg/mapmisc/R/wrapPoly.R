wrapPoly = function(x, crs){
	
	if(is.null(attributes(crs)$crop)) {
		attributes(crs)$crop = llCropBox(crs)
	}
	
	if(requireNamespace('rgeos', quietly=TRUE) & 
			requireNamespace('rgdal', quietly=TRUE)) {	
		toCropX = spTransform(attributes(crs)$crop, crs(x))
		xCrop = rgeos::gDifference(x, toCropX, byid=TRUE)
		xCropData = x@data[match(
						gsub(" [[:digit:]]+$","", names(xCrop)),
						rownames(x@data)
				),]
		rownames(xCropData) = names(xCrop)
		
		xCrop = SpatialPolygonsDataFrame(
				xCrop,
				data=xCropData
		)
		
		xTcrop = spTransform(xCrop, crs)
		
	} else {
		xTcrop = NULL
	}
	
	xTcrop
}

llCropBox = function(crs, 
		crsSphere = mapmisc::crsLL, 
		keepInner=FALSE, N=200) {

	polarLimit = 90
	extentLL = extent(-180,180,-polarLimit,polarLimit)
	if(length(grep("proj=moll", as.character(crs)))){
		
		buffer = 360/N
		
		lon0 = as.numeric(gsub(
						"lon_wrap=", "", 
						regmatches(
								as.character(crs),
								regexpr("lon_wrap=([[:digit:]]|\\.)+ ", 
										as.character(crs))
			)
		))

		lonSplit = Arg(exp(2i*pi*lon0/360 + 1i*pi))*360/(2*pi) 

		lonSeq = exp(seq(0, log(polarLimit), len=N))
		lonSeq = sort(unique(c(lonSeq, -lonSeq)))
		lonMat = rbind(
				cbind(pmax(-180,lonSplit-buffer), lonSeq),
				cbind(pmin(lonSplit+buffer),rev(lonSeq))
				)
		toCropPoly = SpatialPolygons(list(
						Polygons(list(
										Polygon(lonMat, hole=FALSE)
										), 1)
						),
						proj4string = crsSphere)  	
		rasterSphere = rasterize(toCropPoly,
				raster(extentLL, ncol=N, nrow=N, crs=crsSphere),
				field=0
				)
		rasterLL = projectRaster(rasterSphere, crs=mapmisc::crsLL, method='ngb')		
		values(rasterLL)[is.na(values(rasterLL))] = 1
	} else {
	
	Ny = N
	
	rasterLL = raster(
  		extentLL,
			ncol=Ny, nrow=Ny, crs=crsSphere
	)
	rasterT = projectExtent(rasterLL, crs)
	
	rasterTsmall = raster::crop(rasterT, extend(extent(rasterT), -2*res(rasterT)))
	values(rasterTsmall) = 1
	rasterT = extend(rasterTsmall, extent(rasterT), value=0)
	
	rasterLL = projectRaster(rasterT, crs=mapmisc::crsLL, method='ngb')
	if(keepInner){
		values(rasterLL)[is.na(values(rasterLL))] = 1
	} else {
		values(rasterLL)[is.na(values(rasterLL))] = 0
	}
	}
	
	borderLL = rasterToPoints(rasterLL, fun=function(x) {x<1}, spatial=TRUE)
	
	if(requireNamespace('rgeos', quietly=TRUE)) {
		crs(borderLL) = NA
		toCropLL = rgeos::gBuffer(borderLL, width=mean(res(rasterLL)*2))
		crs(toCropLL) = mapmisc::crsLL
	} else {
		toCropLL = NULL
	}
	
	if(!is.null(toCropLL)){
		toCropLL = raster::crop(toCropLL, extentLL)
	}

	toCropLL
}
