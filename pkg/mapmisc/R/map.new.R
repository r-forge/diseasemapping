map.new = function(x, rasterLegend=FALSE) {
	

	
	xpoints = t(bbox(extent(x)))

	xpoints = SpatialPoints(xpoints)
	thecrs = try(proj4string(x), silent=TRUE) 
	if(class(thecrs)!="try-error")
		proj4string(xpoints) = CRS(thecrs)

	
	if(rasterLegend) { 
		par(mar=c(0,0,0,5))
	} else {
		par(mar=c(0,0,0,0))
	}
	
	plot(xpoints,pch=NA)

	
	return(invisible())
		
}
