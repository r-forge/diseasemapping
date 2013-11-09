colourScale = function(x, breaks=5, 
		style=c("quantile","equal","unique", "fixed"),
		col="YlOrRd", opacity=1, dec=NULL, firstBreak=NULL, 
		transform=NULL, revCol=FALSE, exclude=NULL, ...) {
	
	style = style[1]
	if(!is.function(col)){		
		colString = col
		if(length(colString)==1){
			col = function(n) brewer.pal(n, colString)[1:n]
		} else {
			col = function(n) colString[1:n]
		}
	}
	eps=0.01
	
	toexclude = x %in% exclude
	toexclude[is.na(toexclude)] = FALSE
	if(any(toexclude))
		x=x[toexclude]	= NA
	

	
	# unique breaks	
	if(is.character(x) | style=="unique" ) {
		thetable = sort(table(x))[1:breaks]
		breaks = names(thetable)
		colVec = col(length(breaks))
		names(colVec) = breaks
		x = colVec[as.character(x)]		
 
	} else {
		
		if(style != "fixed") {
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
						if(any(x<=0) ) {
							warning("negative or zero x's given with log transform")
						}
						tranform = list(log,exp)
					} else if(transform=="sqrt") {
						if(any(x<=0) ) {
							warning("negative or zero x's given with root transform")
						}
						transform = list(sqrt, function(x) x^2)						
					} else if(transform=="square") {
						if(any(x<=0) ) {
							warning("negative or zero x's given with square transform")
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
				if(is.null(firstBreak))
					firstBreak = min(x, na.rm=TRUE)
				breaks = seq(firstBreak, max(x, na.rm=TRUE),len=breaks)
				firstBreak = NULL
			} else {
				breaks = classIntervals(x, n=breaks, 
						style=style, ...)$brks
			}
			
			
			
			if(!is.null(transform)) {
				breaks = transform[[2]](breaks)
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
		if(revCol) colVec = rev(colVec)
		
		if(any(opacity < 1-eps)) {
			
		}
		
		
		x = as.character(cut(
						x, breaks=breaks,
						labels=colVec,
						include.lowest=TRUE
						))
	}
	
	result = list(col=x,breaks=breaks,legendCol=colVec)
	
		
}
