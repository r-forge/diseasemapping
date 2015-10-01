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

llCropBox = function(crs,N=200) {
	
	Ny = N
	
	rasterLL = raster(extent(-180,180,-90,90), ncol=Ny, nrow=Ny, crs=mapmisc::crsLL)
	rasterT = projectExtent(rasterLL, crs)
	rasterTsmall = raster::crop(rasterT, extend(extent(rasterT), -3*res(rasterT)))
	values(rasterTsmall) = 1
	rasterT = extend(rasterTsmall, extent(rasterT), value=0)
	
	rasterLL = projectRaster(rasterT, rasterLL, method='ngb')
	values(rasterLL)[is.na(values(rasterLL))] = 0
	borderLL = rasterToPoints(rasterLL, fun=function(x) {x<1}, spatial=TRUE)
	crs(borderLL) = NA
	
	if(requireNamespace('rgeos', quietly=TRUE)) {
		toCropLL = rgeos::gBuffer(borderLL, width=mean(res(rasterLL)*1.5))
		crs(toCropLL) = mapmisc::crsLL
		toCropLL = raster::crop(toCropLL, extent(-180,180,-90,90))
	} else {
		toCropLL = NULL
	}
	toCropLL
}
