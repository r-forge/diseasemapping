wrapPoly = function(x, crs){
	
	if(is.null(attributes(crs)$crop)) {
		attributes(crs)$crop = llCropBox(crs)$crop
	}
	
	if(requireNamespace('rgeos', quietly=TRUE) & 
			requireNamespace('rgdal', quietly=TRUE)) {	
		toCropX = spTransform(attributes(crs)$crop, crs(x))
		xCrop = rgeos::gDifference(x, toCropX, byid=TRUE)
		
		if(any(slotNames(x)=='data')) {
		
		xCropData = x@data[match(
						gsub(" [[:digit:]]+$","", names(xCrop)),
						rownames(x@data)
				),]
		rownames(xCropData) = names(xCrop)
		
		xCrop = SpatialPolygonsDataFrame(
				xCrop,
				data=xCropData
		)
		
	}
	
		xTcrop = spTransform(xCrop, crs)
		
	} else {
		xTcrop = NULL
	}
	
	xTcrop
}

llCropBox = function(crs, 
		res=0.5, keepInner=FALSE) {
	
	polarLimit = 90
	extentLL = extent(-180,180,-polarLimit,polarLimit)
	
	
	
	if(length(grep("proj=moll", as.character(crs)))){
		
		projMoll = CRS("+proj=moll +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs84=0,0,0,0,0,0,0 ")
		
		N = round(180/res)
		buffer = 1
		
		lonSeq = exp(seq(0, log(polarLimit-1), len=N))
		lonSeq = sort(unique(c(lonSeq, -lonSeq)))

		lonMat = rbind(
				cbind(-180+1, lonSeq),
				cbind(180-1,rev(lonSeq))
		)
		edgePoints = SpatialPoints(lonMat, proj4string=crsLL)
		edgePointsT = spTransform(edgePoints, projMoll)
		
		toCropPoly = SpatialPolygons(list(
				Polygons(list(
								Polygon(edgePointsT@coords*0.99, hole=FALSE)
						), 1)
		), proj4string = crs)
		rasterT = rasterize(toCropPoly, 
				raster(extent(edgePointsT), nrows=N, ncols=2*N), 
				value=1)
		values(rasterT)[is.na(values(rasterT))]=0
	} else {

		toCropPoly = NULL
		
		rasterLLorig = raster(
  			extentLL,
				res=res, crs=mapmisc::crsLL
		)
		
		rasterTorig = projectExtent(rasterLLorig, crs)
		rasterTorig = disaggregate(rasterTorig, 2)
		
		# put 0's around the border
		rasterTsmall = crop(rasterTorig, extend(extent(rasterTorig), -6*res(rasterTorig)))
		values(rasterTsmall) = 1
		rasterT = extend(rasterTsmall, extend(extent(rasterTsmall), 5*res(rasterTsmall)), value=0)
		
	}
	
	
	rasterLL = projectRaster(
			from=rasterT, 
			crs=mapmisc::crsLL, 
			res = res, method='ngb')
	rasterLL = crop(rasterLL, extentLL)
	if(keepInner){
		values(rasterLL)[is.na(values(rasterLL))] = 1
	} else {
		values(rasterLL)[is.na(values(rasterLL))] = 0
	}
 
	borderLL = rasterToPoints(rasterLL, fun=function(x) {x<1}, spatial=TRUE)
	if(requireNamespace('rgeos', quietly=TRUE)) {
		crs(borderLL) = NA
		toCropLL = rgeos::gBuffer(borderLL, width=mean(res*1.5))
		crs(toCropLL) = mapmisc::crsLL
	} else {
		toCropLL = NULL
	}
	
	if(!is.null(toCropLL)){
		toCropLL = raster::crop(toCropLL, extentLL)
	}

	list(crop=toCropLL, poly=toCropPoly)
}
