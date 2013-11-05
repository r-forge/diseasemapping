
scalebar = function(crs, pos="bottomright",scale.cex=1,...) {
	
	if(is.character(crs))
		crs = CRS(crs)
	if(class(crs) != "CRS")
		crs = CRS(proj4string(crs))
	

#	dash = "\u2517\u2501\u2501\u2501\u2537\u2501\u2501\u2501\u251B"
	
	xpoints = t(bbox(extent(par("usr"))))
	
	xpoints = rbind(xpoints, 
			centre = apply(xpoints, 2, mean))
	xpoints = rbind(xpoints, 
			right = c(xpoints["max",1], xpoints["centre",2]),
			dashright = xpoints["centre",]+c(strwidth("W"),0))
	xpoints = SpatialPoints(xpoints, proj4string=crs)

	
	xll = spTransform(xpoints, CRS("+init=epsg:4326"))
	

	up = matrix(coordinates(xll)["centre",]+c(0,0.1),ncol=2,
			dimnames=list("up",NULL))
	
	
	xll=rbind(xll, SpatialPoints(up, 
					proj4string=CRS(proj4string(xll))))
	
	
	xpoints2 = spTransform(xll[c("up","centre")], CRS(proj4string(xpoints)))
	thediff=apply(coordinates(xpoints2), 2,diff)
	north=atan(thediff[1]/thediff[2])
	

	dashdist = spDists(xll[c("centre","dashright"),], 
			longlat=TRUE)[1,2]*1000 
	bardist = 	dashdist*6*scale.cex
		
	theb = log10(bardist)
	candidates = 10^c(floor(theb), ceiling(theb))
	candidates = c(candidates[1]*c(1,2,5), candidates[2])
	segdist=candidates[order(abs(candidates - bardist))[1]]
	
	segscale = segdist / dashdist
	


	
	if(segdist >900) {
		lunits="km"
		segdist = segdist / 1000
	} else {
		lunits="m"
	}

	eps = 0.175

	theN = c(0, 0+1i, eps+1i, (1-eps)+1.5*eps*1i,
			(1-eps)+1i, 1+1i, 1+0i,
			(1-eps)+0i, eps+(1-1.5*eps)*1i,
			eps) + 0.5+0.5*1i
	theHat = c(-0.25+1.1i, 0.5+1.75i, 1.25+1.1i)
	theHat = c(theHat, rev(theHat) + 2*eps*1i)+ 0.5+0.5*1i
	theN =  strwidth("N")*theN
	theHat =  strwidth("N") * theHat
	
	theN = theN * exp(1i*north)
	theHat = theHat * exp(1i*north)
	
	
	
	thelabel = paste(segdist, lunits,sep="")
	
	forLegend = list(...)
	if(is.null(forLegend$col)) {
		forLegend$col="black"
	} else {
		forLegend$col= forLegend$col[1]
	}
	if(is.null(forLegend$cex)) {
		forLegend$cex = 1
	}
	if(is.null(forLegend$pt.cex)) {
		forLegend$pt.cex = forLegend$cex
	}
	if(is.null(forLegend$text.col)) {
		forLegend$text.col = forLegend$col
	}
	if(is.null(forLegend$title.adj)) {
		forLegend$title.adj = 0.1
	}
	
	forLegend$x = pos
	forLegend$lty = 1
	forLegend$pch = NA
	forLegend$seg.len = segscale
	forLegend$title=thelabel
	forLegend$legend = "       "
	forLegend$lwd=3
		
	
	thelegend = do.call(legend, forLegend)
			
	thecentre =  thelegend$text$x + 1i*thelegend$text$y
	 polygon(forLegend$pt.cex*theN +thecentre, 
			 col=forLegend$text.col,border=NA)
	 polygon(forLegend$pt.cex*theHat + thecentre, 
			 col=forLegend$text.col,border=NA)
	 

	 
	 return(invisible())	
}