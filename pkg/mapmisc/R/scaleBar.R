
scaleBar = function(crs, pos="bottomright",
    scale.cex=1,outer=TRUE,...) {

	if(is.character(crs))
		crs = CRS(crs)
	if(class(crs) != "CRS")
		crs = CRS(proj4string(crs))
	

#	dash = "\u2517\u2501\u2501\u2501\u2537\u2501\u2501\u2501\u251B"
	
	oldcex = par("cex")
	forLegend = list(...)
	if(length(forLegend$cex)){
		par(cex=forLegend$cex)
	} 
	forLegend$cex=1

	xpoints = t(bbox(extent(par("usr"))))
	
	dashTemplate = " 2000 km "
	xcentre = apply(xpoints, 2, mean)
	xpoints = rbind(centre=xcentre, 
			dashright = xcentre + c(strwidth(dashTemplate),0)
	)

	
	xpoints = SpatialPoints(xpoints, proj4string=crs)

	if(requireNamespace('rgdal', quietly=TRUE)) {	
		xll = spTransform(xpoints, crsLL)
	} else {
		xll= xpoints
		if(!length(grep("longlat", projection(xpoints))))
			warning('rgdal not intalled, assuming the plot is long-lat')
	}
	

	up = matrix(coordinates(xll)["centre",]+c(0,0.1),ncol=2,
			dimnames=list("up",NULL))
	
	
	xll=rbind(xll, SpatialPoints(up, 
					proj4string=CRS(proj4string(xll))))
	
	
	

	dashdist = spDists(xll[c("centre","dashright"),], 
			longlat=TRUE)[1,2]*1000 
	bardist = 	dashdist*scale.cex
		
	theb = log10(bardist)
	candidates = 10^c(floor(theb), ceiling(theb))
	candidates = c(candidates[1]*c(1,2,5), candidates[2])
	segdist=candidates[order(abs(candidates - bardist))[1]]
	
	segscale = ( strwidth(dashTemplate)/par("cxy")[1] ) *
			segdist / dashdist
	
	
	if(segdist >900) {
		lunits="km"
		segdist = segdist / 1000
	} else {
		lunits="m"
	}

	eps = 0.175
	
	dimIn = par("pin")
	dimUser = par("usr")
	dimUser = c(dimUser[2]-dimUser[1], dimUser[4]-dimUser[3])
	InPerUnit = dimIn/dimUser
	

	theN = c(0, 0+1i, eps+1i, (1-eps)+1.5*eps*1i,
			(1-eps)+1i, 1+1i, 1+0i,
			(1-eps)+0i, eps+(1-1.5*eps)*1i,
			eps) - 0.5 - 0.5*1i
	theHat = c(-0.25+1i, 0.5+1.6i, 1.25+1i)
	theHat = c(theHat, rev(theHat) + 1.5*eps*1i)- 0.5-0.5*1i
	theN =  strwidth("N")*Re(theN) + 
			1i*strwidth("N")*InPerUnit[1]/InPerUnit[2]*
			Im(theN)
	theHat =  strwidth("N") * Re(theHat)+ 
			1i*strwidth("N")*InPerUnit[1]/InPerUnit[2]*
			Im(theHat)
	
	
	
if(scale.cex>0) {	
	thelabel = paste(segdist, lunits,sep="")
} else {
  thelabel=''
}

	
	defaults = list(col='black', 
			xjust=0.7, bg="white",
			x=pos, text.width=strwidth("I"), pt.cex=1)
	
	if(outer){
		par(xpd=TRUE)
		fromEdge = matrix(par("plt"), 2, 2, 
				dimnames=list(c("min","max"), c("x","y")))
		propIn = apply(fromEdge, 2, diff)
		if(is.character(pos)) {
			inset = c(0,0)
			if(length(grep("^bottom", pos))){
				inset[2] = -fromEdge["min","y"]					
			} else if(length(grep("^top", pos))){
				inset[2] = fromEdge["max","y"]-1					
			}
			if(length(grep("left$", pos))){
				inset[1] = -fromEdge["min","x"]					
			} else if(length(grep("right$", pos))){
				inset[1] = fromEdge["max","x"]-1					
			}
			inset = inset/propIn
			defaults$inset = inset+ 0.01
		}
	} else {
		defaults$inset = 0.01
	}

	for(D in names(defaults)) {
		if(is.null(forLegend[[D]]))
			forLegend[[D]] = defaults[[D]]			
	}

	onem = par("cxy")[1]
	defaults = list(text.col = forLegend$col,
			title.adj = 
					0.5*(segscale*onem/
						(segscale*onem + 4*onem) 
						)* 
					(1-forLegend$text.width/(segscale*onem))
						)
	for(D in names(defaults)) {
		if(is.null(forLegend[[D]]))
			forLegend[[D]] = defaults[[D]]			
	}
	
	forLegend$lty = as.integer(segscale>0)
	forLegend$pch = NA
	forLegend$seg.len = segscale
	forLegend$title=thelabel
	forLegend$legend = NA
	forLegend$lwd=3

	
	
	if(forLegend$seg.len*onem < strwidth(forLegend$title)) {
		forLegend$title=NA
	}

		
	thelegend = do.call(graphics::legend, forLegend)
			
	if(is.na(forLegend$title))
		text(thelegend$text$x - (2/3)*strwidth("m")*forLegend$seg.len,
				thelegend$rect$top , 
				label=thelabel, pos=1, cex=0.75, offset=1.25)
		
  # if there's no scale bar or box and pos is numeric,
# put the N at the point spcified
  if(scale.cex<=0 & !nchar(forLegend$title) & all(forLegend$bty=='n') & is.numeric(forLegend$x)) {	
    thecentre =  c(forLegend$x,forLegend$y)[1:2]
  } else {
    thecentre =  c(thelegend$text$x,thelegend$text$y)
  }
  
  xpoints = SpatialPoints(t(thecentre),
      proj4string=crs)
  
  if(requireNamespace('rgdal', quietly=TRUE)) {	
    xll = spTransform(xpoints, crsLL)
  } else {
    xll= xpoints
    if(!length(grep("longlat", projection(xpoints))))
      warning('rgdal not intalled, assuming the plot is long-lat')
  }
  xll = rbind(xll,
      SpatialPoints(
          xll@coords+c(0,1),
      proj4string=crs(xll))
  )
  
  if(requireNamespace('rgdal', quietly=TRUE)) {	
    xpoints2 = spTransform(xll, crs)
  } else{
    xpoints2 = xll
  }
  thediff=apply(coordinates(xpoints2), 2,diff)
  north=atan(thediff[1]/thediff[2])+pi*(thediff[2]<0)
  
  theN = theN * exp(-1i*north)
  theHat = theHat * exp(-1i*north)
  
  thecentre = thecentre[1] + 1i*thecentre[2]
	 polygon(forLegend$pt.cex*theN +thecentre, 
			 col=forLegend$text.col,border=NA)
	 polygon(forLegend$pt.cex*theHat + thecentre, 
			 col=forLegend$text.col,border=NA)
	 
	par(cex=oldcex)
	 
	return(invisible(list(
              out=thelegend, 
              call=forLegend, 
              centre=c(Re(thecentre),Im(thecentre))
  )))
}