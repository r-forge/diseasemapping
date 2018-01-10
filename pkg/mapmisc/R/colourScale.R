roundForBreaks = function(breaks, dec){
	if(is.null(dec)) return(breaks)
		dec = 10^dec
	breaks = breaks * dec
	breaks = c(floor(breaks[1]), 
		round(breaks[seq(2,by=1,len=length(breaks)-2)]),
		ceiling(breaks[length(breaks)]))
	breaks / dec
}



signifForBreaks = function(breaks, digits=6){
	if(is.null(digits)) return(breaks)

	breaksOrig = breaks
	breaks = signif(breaks, floor(digits))

	if(digits > floor(digits)) {
		breaksBeforeExtra = breaks
		digitsDec = digits - floor(digits)
		toRound = breaksOrig - breaks
		toMult = pmax(1,10^ceiling(log10(abs(toRound))))
		breaksExtra = toRound /toMult
		breaksExtra2 = digitsDec * round(breaksExtra/digitsDec)
		breaks = breaks + breaksExtra2*toMult

		if(breaks[1] > breaksOrig[1])
			breaks[1] = breaksBeforeExtra[1]
		if(breaks[length(breaks)] < breaksOrig[length(breaks)])
			breaks[length(breaks)] = breaksBeforeExtra[length(breaks)]
	}


	toAdd = 10^floor(log10(abs(breaks)))
	if(breaks[1] > breaksOrig[1]) {
		breaks[1] = breaks[1] - toAdd[1]
	}
	if(breaks[length(breaks)] < breaksOrig[length(breaks)]) {
		toAdd = sort(unique(toAdd))
		toAdd = toAdd[pmax(1,length(toAdd)-1)]
		breaks[length(breaks)] = ceiling(
			breaksOrig[length(breaks)]/toAdd) * toAdd
	}
	breaks
}

radarCol <- c("#FFFFFF", "#99CCFF", "#0099FF", "#00FF66", "#00CC00", "#009900", 
	"#006600", "#FFFF33", "#FFCC00", "#FF9900", "#FF6600", "#FF0000", 
	"#FF0299", "#9933CC", "#660099")

colourScale = function(x=NULL, breaks=5, 
	style=c("quantile","equal","unique", "fixed"),
	col="YlOrRd", opacity=1, 
	dec=NULL, digits=6,
	firstBreak=NULL, 
	transform=NULL, revCol=FALSE, exclude=NULL, 
	labels=NULL,...) {


	UseMethod("colourScale")	

}

colorScale = function(...) {
	colourScale(...)
}

colourScale.character =   function(x=NULL, breaks=5, 
	style=c("quantile","equal","unique", "fixed"),
	col="YlOrRd", opacity=1, dec=NULL, digits=6,
	firstBreak=NULL, 
	transform=NULL, revCol=FALSE, exclude=NULL, labels=NULL,...) {

	colourScale(factor(x), breaks=breaks, 
		style=style,
		col=col, opacity=opacity, dec=dec, 
		digits = digits,
		firstBreak=firstBreak, 
		transform=transform, revCol=revCol, exclude=exclude, labels=labels,...)
}

colourScale.factor = function(x=NULL, breaks=5, 
	style=c("quantile","equal","unique", "fixed"),
	col="YlOrRd", opacity=1, dec=NULL, digits=6, firstBreak=NULL, 
	transform=NULL, revCol=FALSE, exclude=NULL, labels=NULL,...) {

	res=colourScale(x=as.integer(x), 
		breaks=breaks,
		style="unique",
		col=col, opacity=opacity, dec=dec, 
digits = digits,
		firstBreak=firstBreak, 
		transform=transform, revCol=revCol, exclude=exclude, labels=as.character(x),
		...)

	res
}

colourScale.Raster = function(x=NULL, breaks=5, 
	style=c("quantile","equal","unique", "fixed"),
	col="YlOrRd", opacity=1, dec=NULL, digits=6,
	firstBreak=NULL, 
	transform=NULL, revCol=FALSE, exclude=NULL, labels=NULL,...) {

	style = style[1]
	weights=NULL

	NforSample = 1e+05

	if(is.character(col) ){

		if(identical(col[1],'radar')){
			col = radarCol
		}
	}


	if(style == "equal") {
		if(length(exclude)) {
			x = unique(x)
		} else {
			x = try(range(c(minValue(x), maxValue(x))), silent=TRUE)
			if(class(x)=='try-error')
				x = range(quantile(sampleRegular(x, size=NforSample)))
		}
	} else if(style=='fixed') {
		x = NULL
	} else if(style=='unique') {

	   # if labels is missing, take labels from the raster
		if(is.null(labels)){
			labels = levels(x)[[1]]
		}

    # if labels is a data frame, use it
		if(is.data.frame(labels)) {
			levelsx = labels
				  # ID variable is first column if it's missing
			if(!any(colnames(levelsx)=='ID')){
				colnames(levelsx)[1] = 'ID'
			}

				  # if a factor or characters, convert to numeric
			if(!is.numeric(levelsx$ID)) {
				levelsx$ID = as.numeric(gsub(
					"^[[:alpha:]]+", "", as.character(levelsx$ID)
					))
			}


				  # different spellings of 'label'
			if(!any(colnames(levelsx)=='label')){
				labelCol = grep("^label(s?)([[:space:]]?)|^NAME$", 
					names(levelsx), ignore.case=TRUE)
				if(length(labelCol)){
					labelCol = labelCol[1]
				} else {
					labelCol = 2
				}
				colnames(levelsx)[labelCol] = 'label'
			}
				  # convert factor or numbers to character
			levelsx$label = as.character(levelsx$label)

			rgbCols = mapply(grep, 
				pattern=paste('^', c("red",'green','blue'), '$', sep=''),
				MoreArgs = list(x=names(levelsx), ignore.case=TRUE)
				)
				  # if all three (rgb) are found
			if(is.numeric(rgbCols)){
				col = levelsx[,rgbCols, drop=FALSE]
				breaks = nrow(col)
			}

			colCol = grep("^col$", colnames(levelsx), ignore.case=TRUE)
			if(length(colCol)) {
				col = levelsx[,colCol]
				breaks = length(col)
			}
		} else if( # labels not a data frame
			length(labels) == length(breaks)
			){
			levelsx = data.frame(
				ID=breaks,
				label=as.character(labels)
				)
		} else { # different number of labels and breaks
		levelsx = data.frame(
			ID=sort(unique(x))
			)
		if(ncol(levelsx)==length(labels)) {
			levelsx$label = as.character(labels)
		} else {
			warning('labels must be same length as either breaks or unique(x)')
			levelsx$label = as.character(levelsx$ID)
		}
	} # end different numbers of levels and breaks


	if(ncell(x)<1e+06) {
		x = freq(x, useNA='no')
		weights = x[,2]
		x=x[,1]
	} else {
		weights = table(
			sampleRegular(x, NforSample)
			)
		x = as.numeric(names(weights))
		if(!is.null(levelsx)) {
				    # check for values which are missing
			toAdd = setdiff(levelsx$ID, x)
			toAddWeights = rep(0, length(toAdd))
			names(toAddWeights) = as.character(toAdd)
			x = c(x, toAdd)
			weights = c(weights, toAddWeights)
		}
	}

	if(is.null(levelsx)){
		levelsx = data.frame(
			ID=sort(x))
	}
	notInLevels = which(! x %in% levelsx$ID)
	if(length(notInLevels)){
      # add more values to ID
		toAdd = matrix(NA,
			length(notInLevels),ncol(levelsx), 
			dimnames=list(NULL, colnames(levelsx)))
		toAdd[,1] = x[notInLevels]
		levelsx = rbind(levelsx, toAdd)
		levelsx = levelsx[order(levelsx$ID),]
	} # end add more ID in levels
	if(is.null(levelsx$label)) {
		if(is.vector(labels)){
			if(length(labels)==nrow(levelsx))
				levelsx$label = labels
		} # end labels are vector
	}
	if(is.null(levelsx$label)) {
		levelsx$label = as.character(levelsx$ID)
	}

	levelsx$freq = weights[match(
		levelsx$ID,
		x
		)]
    # if more breaks than ID's have been requested
    # set breaks to all ID's
	if(length(breaks)==1 & all(breaks > nrow(levelsx)))
		breaks = levelsx$ID

	weights = levelsx$freq
	x=levelsx$ID
	labels = levelsx$label

} else { # not unique or equal or fixed, take a sample
Nsample = 20000
xVec= c()
Diter = 0
while(length(xVec)<Nsample & Diter < 5){
	xVec = c(xVec,
		na.omit(sampleRegular(x, min(c(ncell(x),Nsample))))
		)
	Diter = Diter + 1
}
x = c(xVec, maxValue(x), minValue(x))
} # end if style== stuff

res=colourScale(x, breaks, 
	style,
	col=col, opacity=opacity, dec=dec, digits = digits,
	firstBreak=firstBreak, 
	transform=transform, revCol=revCol, exclude=exclude, labels=labels,
	weights=weights,...)
res[!names(res)%in% 'plot']
}



colourScale.numeric = function(x=NULL, breaks=5, 
	style=c("quantile","equal","unique", "fixed"),
	col="YlOrRd", opacity=1, dec=NULL, digits=6,
	firstBreak=NULL, 
	transform=NULL, revCol=FALSE, exclude=NULL, labels=NULL, 
	weights=NULL,...) {

	xOrig = x
	style = style[1]
	eps=0.01

	 # radar colours
	if(is.character(col)){
		if(identical(col[1],'radar')){
			col = radarCol
		}
	}

	 # convert rgb colours to character string
	if(is.matrix(col)| is.data.frame(col)) {
		redCol = grep("^[[:space:]]*r(ed)?[[:space:]]*$", colnames(col), value=TRUE,ignore.case=TRUE)
		blueCol = grep("^[[:space:]]*b(lue)?[[:space:]]*$", colnames(col), value=TRUE,ignore.case=TRUE)
		greenCol = grep("^[[:space:]]*g(reen)?[[:space:]]*$", colnames(col), value=TRUE,ignore.case=TRUE)
		if(!all(c(length(redCol),length(greenCol), length(blueCol))==1)) {
			warning("col is a matrix but it's not clear which columns are red, green and blue")
		}		
		if(any(col[,c(redCol,greenCol,blueCol)] > 2 )) {
			theMax = 255
			col = round(as.matrix(col[,c(redCol,greenCol,blueCol)]))
		} else{
			theMax = 1
		}
		col = grDevices::rgb(col[,redCol], col[,greenCol], col[,blueCol], maxColorValue=theMax)
	}


	if(!is.function(col)){		
		colName = col
		colString = NULL
		if(length(col)==1){
			if(requireNamespace('RColorBrewer',quietly=TRUE)) {
				col = function(n) RColorBrewer::brewer.pal(n, colName)[1:n]
			} else {
				col = function(n) heat.colors(n)
			}
		} else {
			colString = col
			colName = NULL
			col = function(n) colString[1:n]
		}
	} else {
		colString = NULL
		colName = NULL
	}



	if(style=='unique'){
		  # unique breaks	

		if(length(weights)!=length(x)) {
			thetable = as.data.frame(table(ID=x))
			thetable$ID = as.numeric(as.character(thetable$ID))
		} else {
			thetable = tapply(weights, x,sum, na.rm=TRUE)
			thetable = data.frame(ID=as.numeric(names(thetable)),
				Freq=thetable)
		}

		  # put labels in table
		if(length(labels) == length(breaks)) {
			   # labels are assigned to breaks
			thetable$label = labels[match(thetable$ID,breaks)]
		} else if (length(labels) == length(x)) {
			   # labels algined with x
			thetable$label = labels[match(thetable$ID, x)]
		} else if(length(labels) == nrow(thetable)) {
			   # labels are assigned to x
			thetable$label = labels
		} else {
			   # don't use labels
			thetable$label = as.character(thetable$ID)
		}

		if(length(breaks) > 1) {

		    # breaks not in the table
			notInX = which(! breaks %in% thetable$ID)
			if(length(notInX))
				thetable = rbind(thetable,
					data.frame(ID=breaks[notInX],
						Freq=rep(0,length(notInX)),
						label = as.character(breaks[notInX]))
					)
		}

		  # if colours is a vector see if it can be put into table
		if (length(breaks) == length(colString)) {
			thetable$col = colString[match(thetable$ID,breaks)]
		} else if(nrow(thetable) == length(colString)) {
			thetable$col = colString
		} 

		shouldExclude = which(thetable$ID %in% exclude)
		if(length(shouldExclude))
			thetable[shouldExclude,'Freq'] = -1


		if(length(breaks)==1){
      # breaks is the maximum number of breaks
			thetable = thetable[order(thetable$Freq, decreasing=TRUE),]

			ncol = min(
				breaks,sum(thetable$Freq>=0) 
				) 
			breaks = thetable[which(thetable$Freq>0)[1:ncol],	'ID']
		}  else {
			breaks = breaks[! ( breaks %in% exclude) ]
		}

		  # assign colours if they're not yet in the table
		if(!length(thetable$col)) {
			thetable$col = NA
			colString = col(length(breaks))
			thetable[match( breaks, thetable$ID),'col'] = colString
		} else {
			   # remove colours from excluded categories
			thetable[thetable$Freq<0, 'col'] = NA
		}

		thetable = thetable[order(thetable$ID),]

		colVec = thetable$col
		names(colVec) =as.character(thetable$ID)
		breaks = thetable$ID
		breaks = c(breaks[1]-1/2, breaks+c(diff(breaks)/2,1/2))
	} else { # not unique breaks


	if(length(exclude) & length(x)) {
		toexclude = which(x %in% exclude)
		x[toexclude]	= NA
	}


	thetable=NULL
	if(!is.null(transform)) {
		if(is.numeric(transform)) {
				    # assume it's a box=cox parameter
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
		} # end transform numeric 
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
		} # end transform character 	
	} else {# end transform not null
	transform = list(function(x) x)[c(1,1)]		
}
xOrig = x
x = transform[[1]](xOrig)

if(style=="quantile"){
	breaks = quantile(x, prob=seq(0,1,len=breaks[1]), na.rm=TRUE)
} else if(style=="equal"){
	startHere = min(x, na.rm=TRUE)
	if(is.null(transform) & !is.null(firstBreak) )	{
		startHere = firstBreak
		firstBreak = NULL
	}					
	breaks = seq(startHere, max(x, na.rm=TRUE),len=breaks[1])
} else if(style!='fixed') {
		    # style is passed to classInt
	if (requireNamespace("classInt", quietly = TRUE)) { 
		if(length(x)>2100)
			x = sample(length(x), 2000)
		breaks = classInt::classIntervals(x, n=breaks, 
			style=style, ...)$brks
	} else {
		warning("Install the classInt package to use style=", style)
	}
} # end classint

breaks = transform[[2]](breaks)

# round
if(!is.null(dec)) {
	breaks = roundForBreaks(breaks, dec)
} 
if(!is.null(digits)) {
	breaks = signifForBreaks(breaks, digits)
} # end rounding

if(!is.null(firstBreak))
	breaks[1] = firstBreak

breaks = sort(unique(breaks))			

colVec = col(length(breaks)-1)		
if(revCol) colVec = rev(colVec)
} # end style not unique


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

result = list(col=colVec, breaks=breaks, colOpacity=colForPlot)
if(style=="unique") {
	thetable$colOpacity = colForPlot

	result$colourtable = rep(NA, max(thetable$ID)+1)
	result$colourtable[1+thetable$ID] = thetable$colOpacity
	result$colortable = result$colourtable

	if(length(xOrig))
		result$plot = thetable[match(xOrig, thetable$ID),'colOpacity']
	result$levels = thetable
	if(length(thetable$label))
		result$legend = as.character(thetable$label)
} else if (length(xOrig)){		
	result$plot = as.character(cut(
		xOrig, 
		breaks=breaks,
		labels=colForPlot,
		include.lowest=TRUE
		))
	if(length(labels)==length(result$col))
		result$legend = labels
}

result

}

colourScale.NULL = colourScale.numeric



colourScale.logical = function(x=NULL, breaks=5, 
	style=c("quantile","equal","unique", "fixed"),
	col="YlOrRd", opacity=1, dec=NULL, digits=6, firstBreak=NULL, 
	transform=NULL, revCol=FALSE, exclude=NULL, labels=NULL, ...) {

	colourScale(as.numeric(x), breaks=breaks, 
		style=style,
		col=col, opacity=opacity, dec=dec, firstBreak=firstBreak, 
		transform=transform, revCol=revCol, exclude=exclude, labels=labels,...)

}
