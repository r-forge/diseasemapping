
scaleBar = function(crs, 
    pos="bottomright",
    cex=1,
    pt.cex = 1.1*cex,
    seg.len=5*cex,
    title.cex=cex,
    outer=TRUE, ...) {

  forLegend = list(cex=cex,pt.cex=pt.cex,seg.len=seg.len,...)
  
  if(length(forLegend$scale.cex)){
    warning('use seg.len instead of scale.cex')
    forLegend$seg.len = forLegend$scale.cex
    forLegend$scale.cex=NULL
  }
  
  
  # if text.width=0, no north arrow
  if(pt.cex==0) {
    forLegend$pt.cex = 0.1*forLegend$cex
    noArrow=TRUE
  } else {
    noArrow=FALSE
  }
  if(seg.len==0){
    forLegend$seg.len = 0.1
    noScale = TRUE
  } else {
    noScale = FALSE
  }

  if(!length(forLegend$text.width)){
    forLegend$text.width = 2*strwidth('m', cex=forLegend$pt.cex)
  }
  

	

#	dash = "\u2517\u2501\u2501\u2501\u2537\u2501\u2501\u2501\u251B"

  # we'll target a label this wide
  dashTargetWidth = strwidth('m', cex=forLegend$cex)*forLegend$seg.len
  
  
	xpoints = matrix(par('usr'), ncol=2)
	xcentre = apply(xpoints, 2, mean)

	xpoints = vect(
			rbind(
				centre=xcentre, 
				dashright = xcentre + c(dashTargetWidth,0)
			), crs=crs(crs))


		xll = terra::project(xpoints, crsLL)

	
# how long (in m) is our target dash
	dashTargetDist = as.matrix(terra::distance(xll))[1,2]
	desiredLen = strwidth("m")*forLegend$seg.len 
	
	theb = log10(dashTargetDist)
	candidates = 10^c(floor(theb), ceiling(theb))
	candidates = c(candidates[1]*c(1,2,5), candidates[2])
	dashRoundedDist = candidates[
			order(abs(candidates - dashTargetDist))[1]
	]
	dashRoundedWidth = dashTargetWidth * dashRoundedDist /
			dashTargetDist
	
	# seg.len gets multipled by 
	# cex * par("cex") * xinch(par('cin')[1L], warn.log = FALSE)
	# in graphics::legend
	forLegend$seg.len = dashRoundedWidth / (
				forLegend$cex * par('cex') * graphics::xinch(par('cin')[1L], warn.log = FALSE)
				)
		
	if(dashRoundedDist >1100) {
		lunits="km"
		segdistPrint = dashRoundedDist / 1000
	} else {
		lunits="m"
    segdistPrint = dashRoundedDist
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
  
  Nwidth = strwidth('N', cex=forLegend$pt.cex)
	theHat = c(theHat, rev(theHat) + 1.5*eps*1i)- 0.5-0.5*1i
	theN =  Nwidth*Re(theN) + 
			1i*Nwidth*InPerUnit[1]/InPerUnit[2]*
			Im(theN)
	theHat =  Nwidth* Re(theHat)+ 
			1i*Nwidth*InPerUnit[1]/InPerUnit[2]*
			Im(theHat)
	
	
	
if(!noScale) {	
	thelabel = paste(segdistPrint, lunits,sep="")
} else {
  thelabel=''
}

	
	defaults = list(col='black', 
			xjust=0.7, 
			x=pos)
	
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
  if(!length(forLegend$text.col))
    forLegend$text.col = forLegend$col

	onem = par("cxy")[1]
	
	forLegend$lty = as.integer(!noScale)
	forLegend$title=' '
	forLegend$legend = ' '
	forLegend$lwd=3

	if(par('bg')=='transparent' & !any(c('bg','bty') %in% names(forLegend))) forLegend$bty = 'n'
	
	thelegend = do.call(graphics::legend, forLegend)
  
  thelegend$title = thelabel
  
  thelegend$textxy =  
       mean(c(thelegend$rect$left, thelegend$text$x)) + 
          1i*thelegend$text$y
			
      
      
	if(forLegend$lty)
		graphics::text(Re(thelegend$textxy),
        Im(thelegend$textxy),
				label=thelegend$title,  
        pos=3, cex=title.cex,
        offset=0.5,
        col=forLegend$text.col)
		
  # if there's no scale bar or box and pos is numeric,
# put the N at the point spcified
  if(noScale & !nchar(forLegend$title) & all(forLegend$bty=='n') & is.numeric(forLegend$x)) {	
    thecentre =  c(forLegend$x,forLegend$y)[1:2]
  } else {
    thecentre =  c(thelegend$text$x + forLegend$text.width/2, thelegend$text$y)
  }
  
	xpoints = vect(
			matrix(thecentre, ncol=2), crs=crs(crs))





				up = vect(
						matrix(terra::crds(xpoints)+c(0,0.1),
								ncol=2,
								dimnames=list("up",NULL)),
						crs = crs(xll)
				)
				
		xll = terra::project(xpoints, crsLL)
				

  xll = rbind(xll,
      vect(
          terra::crds(xll)+c(0,1),
      crs=crsLL)
  )
  
   xpoints2 = terra::project(xll, crs)

  thediff=apply(terra::crds(xpoints2), 2,diff)
  north=atan(thediff[1]/thediff[2])+pi*(thediff[2]<0)
  
  theN = theN * exp(-1i*north)
  theHat = theHat * exp(-1i*north)
  
  thecentre = thecentre[1] + 1i*thecentre[2]
  
  if(!noArrow){
	 polygon(forLegend$pt.cex*theN +thecentre, 
			 col=forLegend$text.col,border=NA)
	 polygon(forLegend$pt.cex*theHat + thecentre, 
			 col=forLegend$text.col,border=NA)
  }
 
 
	return(invisible(list(
              out=thelegend, 
              call=c(forLegend, list(label=thelabel)), 
              centre=c(Re(thecentre),Im(thecentre))
  )))
}