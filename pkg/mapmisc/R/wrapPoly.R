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
						gsub(" border$","", names(xCrop)),
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
		buffer = res
		
		lonSeq = exp(seq(0, log(polarLimit-1), len=N))
		lonSeq = sort(unique(c(lonSeq, -lonSeq)))

		N=100
		lonMat = cbind(
				rep(180, N),
				90-c(0, exp(seq(0,log(90),len=N-1)))
		)
		edgePoints = SpatialPoints(lonMat, proj4string=crsLL)
		edgePointsT = spTransform(edgePoints, projMoll)
		
		edgeCoords = abs(edgePointsT@coords)
		edgeCoords = approx(edgeCoords[,1], edgeCoords[,2], 
				xout = max(edgeCoords[,1]) * sin(seq(0, pi/2, len=N)),
				rule=2)
		
		edgeCoords = cbind(
				c(edgeCoords$x[-1], rev(edgeCoords$x)[-1], -edgeCoords$x[-1], -rev(edgeCoords$x)[-1]),
				c(edgeCoords$y[-1], -rev(edgeCoords$y)[-1], -edgeCoords$y[-1], rev(edgeCoords$y)[-1])
				)
				
				
		toCropPoly = SpatialPolygons(list(
				Polygons(list(
								Polygon(edgeCoords, hole=FALSE)
						), 1)
		), proj4string = crs)

	# clip some of the extreme points because they can loop around
	edgeCoords = edgeCoords[
			abs(edgeCoords[,1]) < 0.999*max(abs(edgeCoords[,1])) &
					abs(edgeCoords[,2]) < 0.999*max(abs(edgeCoords[,2])), 			
			]


	edgePoints = SpatialLines(list(Lines(
			Line(edgeCoords[1:floor(nrow(edgeCoords)/2),]),
			ID = "border"
	)), proj4string=crs)
		
		edgePointsLL = spTransform(edgePoints, crsLL)
		
		edgePointsLL = raster::crop(edgePointsLL, extent(-180+buffer, 180-buffer, -90, 90))
		
		edgePointsLL = edgePointsLL@lines[[1]]@Lines

		toCropLL = SpatialPolygons(
				list(Polygons(
								lapply(edgePointsLL, 
										function(qq){
							Polygon(
								rbind(
										qq@coords[1,],
										cbind(qq@coords[,1], qq@coords[,2] + buffer),
										cbind(rev(qq@coords[,1]), rev(qq@coords[,2])-buffer),
										qq@coords[1,]), 
								hole=FALSE)
				}), ID="border")),
			proj4string=crsLL)
		
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
		
	}
	

	list(crop=toCropLL, poly=toCropPoly)
}
