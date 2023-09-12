map.new = function(x, legendRight=FALSE, buffer=0,
                   mar=c(0,0,0,0),  ...) {
	
	if(!is.null(attributes(x)$ellipse)) {
		ellipse = x = attributes(x)$ellipse
	} else {
		ellipse = NULL
	}
	
  xpoints = terra::vect(
  	terra::crds(terra::vect(terra::extend(
  	terra::ext(x), buffer))), 
  crs=terra::crs(x))

  

  if(is.logical(legendRight))
    legendRight = c(1,0.8)[1+ legendRight ]
  
  if(!is.na(legendRight))
    mar[4] = (1-legendRight) * par('fin')[1]/par('csi')

	ldots = list(...)
	ldots$mar = mar
	if(!length(ldots$bty)) {
		if( mar[4]>0 ) {
			ldots$bty='l'
		} else {
			ldots$bty='n'			
		}
	}
	if(!length(ldots$xpd)) ldots$xpd=FALSE
	if(!length(ldots$axes)) ldots$axes=FALSE
	if(!length(ldots$xaxs)) ldots$xaxs = 'i'
	if(!length(ldots$yaxs)) ldots$yaxs = 'i'
	if(!length(ldots$xaxt)) ldots$xaxt = 'n'
	if(!length(ldots$yaxt)) ldots$yaxt = 'n'
	
  oldpar = par(no.readonly = TRUE) 

  do.call(par, ldots[setdiff(names(ldots),'axes')]) #(mar=mar, xaxs='i',yaxs='i',xpd=FALSE, ...)
	do.call(terra::plot, c(list(x=xpoints,pch=NA), ldots))#, axes=TRUE, ...)

	if(!is.null(ellipse)) {
		plot(ellipse, add=TRUE, ...)
	}
  
  return(invisible(oldpar))
		
}
