map.new = function(x, legendRight=FALSE) {
	

	
	xpoints = t(bbox(extent(x)))

	xpoints = SpatialPoints(xpoints)
	thecrs = try(proj4string(x), silent=TRUE) 
	if(class(thecrs)!="try-error") {
		if(!is.na(thecrs))
			proj4string(xpoints) = CRS(thecrs)
	}
		
#	oldpar = par()[c('mar','plt','xpd')]
#	oldpar = par()[c('mar')]
	
	par(mar=c(0,0,0,0))
	if(legendRight) { 
		if(!is.logical(legendRight)) {
			bob=legendRight
		} else {
			bob=0.8
		}
		par(mar=c(0,0,0,0),plt=c(0,bob, 0,1),xpd=FALSE)
	} else {
		par(mar=c(0,0,0,0),plt=c(0,1, 0,1))
	}
	
	plot(xpoints,pch=NA)

# 	par(oldpar)
	
	return(invisible())
		
}
