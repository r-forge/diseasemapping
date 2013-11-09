map.new = function(x) {
	

	
	xpoints = t(bbox(extent(x)))

	xpoints = SpatialPoints(xpoints)
	thecrs = try(proj4string(x), silent=TRUE) 
	if(class(thecrs)!="try-error")
		proj4string(xpoints) = CRS(thecrs)

	
	par(mar=c(0,0,0,0))
	plot(xpoints,pch=NA)

	
	return(invisible())
		
}
