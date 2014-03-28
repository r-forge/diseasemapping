colourScale = function(x, breaks=5, 
style=c("quantile","equal","unique", "fixed"),
col="YlOrRd", opacity=1, dec=NULL, firstBreak=NULL, 
transform=NULL, revCol=FALSE, exclude=NULL, ...) {

UseMethod("colourScale")	

}

colourScale.factor = function(x, breaks=5, 
		style=c("quantile","equal","unique", "fixed"),
		col="YlOrRd", opacity=1, dec=NULL, firstBreak=NULL, 
		transform=NULL, revCol=FALSE, exclude=NULL, ...) {

	toexclude = which(x %in% exclude)
	if(length(toexclude)){	
		x=x[toexclude] = NA
	}
	xlabels = levels(x)
	names(xlabels) = as.character(seq(1:length(xlabels)))
	x = as.integer(x)
	
res=colourScale(x, breaks=breaks,
			style="unique",
			col, opacity, dec, firstBreak, 
			transform, revCol, exclude, ...) 		

	res$breaks = xlabels[res$breaks]
	res
}
colourScale.Raster = function(x, breaks=5, 
style=c("quantile","equal","unique", "fixed"),
col="YlOrRd", opacity=1, dec=NULL, firstBreak=NULL, 
transform=NULL, revCol=FALSE, exclude=NULL, ...) {

style = style[1]
	if(style %in% c("equal","fixed")) {
		x = c(minValue(x), maxValue(x))
	} else {
		if(ncell(x)<4E+05) {
			x = values(x)
		} else {
			stop("x has too many cells for method", style)
		}
	}
	res=colourScale(x, breaks, 
			style,
			col, opacity, dec, firstBreak, 
			transform, revCol, exclude, ...) 		
	res
}

colourScale.numeric = function(x, breaks=5, 
		style=c("quantile","equal","unique", "fixed"),
		col="YlOrRd", opacity=1, dec=NULL, firstBreak=NULL, 
		transform=NULL, revCol=FALSE, exclude=NULL, ...) {
	
	style = style[1]
	if(!is.function(col)){		
		colString = col
		if(length(colString)==1){
			col = function(n) RColorBrewer::brewer.pal(n, colString)[1:n]
		} else {
			col = function(n) colString[1:n]
		}
	}
	eps=0.01
	
	if(missing(x))
		x=NULL

	
	if(length(exclude) & length(x)) {
	toexclude = x %in% exclude
	toexclude[is.na(toexclude)] = FALSE
	if(any(toexclude,na.rm=TRUE))
		x=x[toexclude]	= NA
	}
	

	
	# unique breaks	
	if(is.character(x) | style=="unique" ) {
		thetable = sort(table(x),decreasing=TRUE)
		thetable = thetable[seq(1, min(c(breaks,length(thetable) )))]
		breaks = names(thetable)
		colVec = col(length(breaks))
		names(colVec) = breaks
 
	} else {
		if(style != "fixed" & length(breaks)==1) {
			if(!is.null(transform)) {
				if(is.numeric(transform)) {
					transform = transform[1]
					if(abs(transform)<eps){
						transform="log"
					} else if(abs(transform-1)<eps) {
						transform=NULL
					} else if(abs(transform-0.5)<eps) {
						transform = "sqrt"
					} else if(abs(transform-2)<eps) {
						transform = "square"
					} else {
						if(any(x<0,na.rm=TRUE) ) {
							warning("negative x's given with Box-Cox parameter")
						}
						boxcox = transform		
						transform = list(
								function(x) {
									(-1)^(boxcox<0)*(x^boxcox - 1)/boxcox
								},
								function(x) {
									((-1)^(boxcox<0)*x*boxcox+1)^(1/boxcox)									
								}
								)						
					}
					
					# assume it's a box=cox parameter
				} 

				
				if(is.character(transform)){
					if(transform=="exp") { 
						tranform = list(exp, log)
					} else if(transform=="log") {
						if(any(x<=0,na.rm=TRUE) ) {
							warning("negative or zero x's given with log transform")
						}
						transform = list(log,exp)
					} else if(transform=="sqrt") {
						if(any(x<0,na.rm=TRUE) ) {
							warning("negative x's given with root transform")
						}
						transform = list(sqrt, function(x) x^2)						
					} else if(transform=="square") {
						if(any(x<0,na.rm=TRUE) ) {
							warning("negative x's given with square transform")
						}
						transform = list(function(x) x^2, sqrt)						
					}
				} 

				xOrig = x
				x = transform[[1]](x)
			}
			
			
			if(style=="quantile"){
				breaks = quantile(x, prob=seq(0,1,len=breaks), na.rm=TRUE)
			} else if(style=="equal"){
				startHere = min(x, na.rm=TRUE)
				if(is.null(transform) & !is.null(firstBreak) )	{
						startHere = firstBreak
						firstBreak = NULL
				}					
				breaks = seq(startHere, max(x, na.rm=TRUE),len=breaks)
			} else {
				breaks = classInt::classIntervals(x, n=breaks, 
						style=style, ...)$brks
			}
		
			
			if(!is.null(transform)) {
				breaks = transform[[2]](breaks)
				if(!is.null(firstBreak))
					breaks[1] = firstBreak
				x = xOrig
			}
			# round
			if(!is.null(dec)) {
				dec = 10^dec
				breaks = breaks * dec
				breaks = c(floor(breaks[1]), 
						round(breaks[seq(2,by=1,len=length(breaks)-2)]),
						ceiling(breaks[length(breaks)]))
				breaks = breaks / dec
			}
			if(!is.null(firstBreak))
				breaks[1] = firstBreak
			
			# unique and sort
			 breaks = sort(unique(breaks))			
			 
		}
		
		colVec = col(length(breaks)-1)		
	}
		

	if(revCol) colVec = rev(colVec)
		
	if(any(opacity < 1-eps)) {
			if(length(opacity)==1){
				opacity = rep(opacity, length(colVec))
			} else if(length(opacity)==2){
				opacity = seq(opacity[1], opacity[2], len=length(colVec))
			} else if(length(opacity)==3){
				opacity = c(opacity[1], 
						seq(opacity[2], opacity[3],
								len=length(colVec)-1))
			} else {
				opacity = opacity[round(
								seq(1,length(opacity), 
										len=length(colVec)))]
			}
			opacity = toupper(as.hexmode(round(opacity*255)))
			hasOpacity = grep("^#[[:xdigit:]]{8}$", colVec)
			colVec[hasOpacity] = substr(colVec, 1, 7)
			
			isRGB = grep("^#[[:xdigit:]]{6}$", colVec)
			colForPlot = colVec
			colForPlot[isRGB] = paste(colForPlot[isRGB], opacity[isRGB],sep="")
			
	} else {
			colForPlot = colVec
	}
		
	if(style=="unique") {
		x = colVec[as.character(x)]
	} else if (length(x)){		
		x = as.character(cut(
						as.numeric(x), 
						breaks=breaks,
						labels=colForPlot,
						include.lowest=TRUE
			))
	}
	
	
	result = list(col=colVec, plot=x,
			breaks=breaks, colOpacity=
					colForPlot)
	
	result
		
}
