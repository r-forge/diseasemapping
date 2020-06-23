
gridlinesWrap = function(crs, 
		easts=seq(-180,180,by=60),
		norths=seq(-90,90,by=30),
		ndiscr=40, plotLines=TRUE, 
		plotLabels = TRUE, ...){
	
	if(any(easts==180)) easts = easts[easts > -180]
	norths = norths[abs(norths) < 90]
	
	pointPos = 2 # 1 for x direction
	
	crsT = crs(crs)
	
	
	glines = c(
			mapply(function(qq, len){
						Lines(
								Line(cbind(qq, 
												seq(-89,89,len=len))),
								ID=paste("E",qq,sep=''))
					}, qq=easts, len=ndiscr), 
			mapply(function(qq, len){
						Lines(
								Line(cbind(
												seq(-179,179,len=len), 
												qq)
								),
								ID=paste("N",qq,sep='')
						)
					}, qq=norths, len=ndiscr)
	)
	
	glines = SpatialLines(glines,
			proj4string = crsLL
	)
	
	
	glinesT = wrapPoly(x=glines, crsT)
	
	ellipseSmall = attributes(crsT)$ellipse
	if(!is.null(ellipseSmall)) {
		ellipseSmall@polygons[[1]]@Polygons[[1]]@coords = 
			0.99*ellipseSmall@polygons[[1]]@Polygons[[1]]@coords 
			
		glinesT = rgeos::gIntersection(
		  spTransform(glinesT,crs(ellipseSmall)), 
		  ellipseSmall, byid=TRUE)
	}
	glinesData=data.frame(
			direction = substr(names(glinesT), 1,1),
			degrees = as.numeric(gsub(
							"^(E|N)|( (buffer|[[:digit:]]+))$", "",
							names(glinesT))
			),
			stringsAsFactors=FALSE
	)
	theNeg = glinesData$degrees < 0
	glinesData[theNeg, 'direction'] = 
			c('W', 'S')[1+(glinesData[theNeg, 'direction']=='N')]
	glinesData$degrees = abs(glinesData$degrees)
	rownames(glinesData) = names(glinesT)
	
	glinesT = SpatialLinesDataFrame(
			glinesT,
			data=glinesData
	)
	
	
	fun1 = function(onedir) {		
		lapply(
				onedir@Lines, 
				function(oneline) {
				  result = oneline@coords[
							which.max(abs(oneline@coords[,pointPos]))
							,]
				  result
				}
		)
	}
	
	legendPoints = lapply(glinesT@lines, 
			function(qq2) {
				t(simplify2array(fun1(qq2)))
			}
	)

	# at most three labels per line
	manyPoints = which(unlist(lapply(legendPoints, nrow)) > 3)
	for(D in manyPoints) {
	  legendPoints[[D]] = legendPoints[[D]][c(1, nrow(legendPoints[[D]])), ]
	}
		
	legendDf = glinesT@data[rep(1:length(glinesT),
					unlist(lapply(legendPoints, nrow))
			),]
	legendDf$label = paste(legendDf$degrees, legendDf$direction)
	
	legendPoints  = SpatialPoints(
			do.call(rbind, legendPoints),
			proj4string=crsT
	)
	rownames(legendDf) = names(legendPoints)
	legendPoints = SpatialPointsDataFrame(
			legendPoints, legendDf
	)
	
	if(plotLines){
		lines(glinesT, ...)
	}		
	if(plotLabels){
		text(legendPoints, labels=legendPoints$label, ...)
	}
	
	invisible(list(
					lines = glinesT,
					points = legendPoints
			)
	)
}

