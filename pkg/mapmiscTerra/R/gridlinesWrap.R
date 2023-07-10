
gridlinesWrap = function(crs, 
		easts=seq(-180,180,by=60),
		norths=seq(-90,90,by=30),
		ndiscr=40, plotLines=TRUE, 
		plotLabels = TRUE, ...){
	
	if(any(easts==180)) easts = easts[easts > -180]
	norths = norths[abs(norths) < 90]
	
	pointPos = 2 # 1 for x direction
	
	crsT = crs(crs)
	
	
	vertSeq = seq(-90,90,len=ndiscr)
	vertLines = expand.grid(y=vertSeq, easts = easts)
	vertLines = vect(as.matrix(vertLines[,c('easts','y')]), type='points', crs=crsLL, 
		atts = data.frame(type='east', loc=vertLines[,'easts']))
	vertLines = terra::split(vertLines,'loc')

	horizSeq = seq(-180,180,len=2*ndiscr)
	horizLines = expand.grid(x=horizSeq, norths = norths)
	horizLines = vect(as.matrix(horizLines[,c('x','norths')]), type='points', crs=crsLL, 
		atts =data.frame(type='east', loc=horizLines[,'norths']))
	horizLines = terra::split(horizLines,'loc')

	glines = vect(lapply(c(horizLines, vertLines), as.lines))
	values(glines) = data.frame(type = rep(c('north','east'), c(length(norths), length(easts))), loc = c(norths, easts))
	glines$neg = glines$loc < 0
	glines$degrees = abs(glines$loc)
	glines$direction = c(north="N",east="E")[glines$type]
	glines[glines$type == 'north' & glines$neg, 'direction'] = 'S'
	glines[glines$type == 'east' & glines$neg, 'direction'] = 'W'
	values(glines)$ID = paste0(glines$degrees, glines$direction)

	glinesT = wrapPoly(x=glines, crsT)

	if(!all(c('ellipse','crop') %in% names(attributes(crsT)$ellipse))) {
		ellipseAndCrop = llCropBox(crsT)
		attributes(crsT)$ellipse = ellipseAndCrop$ellipse
		attributes(crsT)$crop = ellipseAndCrop$crop
	} 

	
	ellipseSmall =  attributes(crsT)$ellipse
	if(!is.null(ellipseSmall)) {
		glinesT = terra::crop(glinesT, ellipseSmall)
	}

	
	legendPoints = terra::centroids(glinesT, inside=TRUE)

	# at most three labels per line
	theTable = table(legendPoints$ID)
	manyPoints = names(theTable[theTable > 3])
	okPoints = which(! legendPoints$ID %in% manyPoints)

	toTrim = data.frame(index = which(legendPoints$ID %in% manyPoints))
	toTrim$ID = legendPoints$ID[toTrim$index]

	okPoints = c(okPoints, unlist(lapply(split(toTrim$index, toTrim$ID), function(xx) xx[seq(from=1, len=3, by=pmax(1,floor(length(xx)/3)))] )))

	legendPoints[okPoints]
	legendPoints$minDist = apply(distance(legendPoints, legendPoints), 2, function(xx) min(xx[xx>0]))


	if(plotLines){
		lines(glinesT, ...)
	}		
	if(plotLabels){
		legendPoints$isClose = legendPoints$minDist < strwidth('XX')
		text(legendPoints[!legendPoints$isClose], labels=legendPoints$ID[!legendPoints$isClose], ...)
	}
	
	invisible(list(
					lines = glinesT,
					points = legendPoints
			)
	)
}

