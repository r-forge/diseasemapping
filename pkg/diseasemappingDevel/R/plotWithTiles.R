plotNamesWithTiles = function(x, attr = 1, brks = NULL, prob = FALSE, 
	Ncol = 5, trans = 50, zoom = 4, xlim = NULL, ylim = NULL,
	filename = NULL, devOff = TRUE, width = 600, height = 600, label = 1
)
{
	library(classInt)
	library(webmaps)

	if(prob == TRUE)
	{
		brks <- c(0,0.2,0.8,0.95,1)
		colours <- c("#00FF00","#FFFF00","#FF6600","#FF0000")
	} else{
		library(RColorBrewer)
		colours <- brewer.pal(Ncol, "YlOrRd")
	}

	if(is.null(brks))
	{
		ci <- classInt::classIntervals(x@data[,attr], Ncol, style = "kmeans")
		brks <- ci$brks * 100
		brks <- c(floor(brks[1]), round(brks[-c(1,length(brks))]), ceiling(brks[length(brks)]))/100
	}

	if(class(x) %in% c("SpatialGridDataFrame","SpatialPixelsDataFrame"))
	{
		library(raster)
		x <- raster(x, layer = attr)
		x <- projectRaster(x, crs = "+proj=longlat +datum=NAD83")
	} else {
		stop("x must be SpatialGridDataFrame or SpatialPixelsDataFrame object")
	}
	
	if(is.null(xlim))
	{
		xlim= bbox(x)[1,]
	}
	if(is.null(ylim))
	{
		ylim= bbox(x)[2,]
	}
	
	bgMap <- getTiles(xlim, ylim, zoom, path = "http://a.www.toolserver.org/tiles/osm-no-labels/", maxTiles = 400)
		
	if(!is.null(filename))
	{ 
		if(!length(grep("\\.png$", filename)))
		{
			png(paste(filename, ".png", sep = ""), width = width, height = height)
		}	
	}
	
	plot(1, type = "n", xlim = xlim,ylim=ylim, xlab=NA,ylab=NA,axes=F)

	image(bgMap, add = TRUE)
	
	brkstemp <- brks
	brkstemp[length(brkstemp)] <- brkstemp[length(brkstemp)] + 0.01
	plot(x, breaks = brkstemp, col = paste(colours, substring(hsv(alpha = as.numeric(trans)/100),8), sep = ""), 
		legend = FALSE, horizontal = FALSE, add = TRUE)

	legend("bottomright", fill = colours, cex = label, 
		legend = c(paste("[", brks[1:(length(brks) - 2)], ",", brks[-c(1,length(brks))], ")", sep = ""), 
			paste("[", brks[length(brks) - 1], ",", brks[length(brks)], "]", sep = "")), bg = "white")
	text(-81.23304, 42.98339, "London", cex = label)
	text(-82.40407, 42.97866, "Sarnia", cex = label)
	text(-82.13981, 42.87978, "Petrolia", cex = label)
	lines(c(-82.21666,-81.97337), c(42.55,42.55), lwd = 4)
	lines(c(-82.21666,-82.21666), c(42.54,42.56), lwd = 2)
	lines(c(-81.97337,-81.97337), c(42.54,42.56), lwd = 2)
	text(-82.095015, 42.52, "20km", cex = label)

	if(!is.null(filename) & devOff)
	{
		dev.off()
	}
}

plotNoNamesWithTiles = function(x, attr = 1, brks = NULL, prob = FALSE, 
	Ncol = 5, trans = 50, zoom = 4, xlim = NULL, ylim = NULL,
	filename = NULL, devOff = TRUE, width = 600, height = 600, label = 1
)
{
	library(classInt)
	library(webmaps)

	if(prob == TRUE)
	{
		brks <- c(0,0.2,0.8,0.95,1)
		colours <- c("#00FF00","#FFFF00","#FF6600","#FF0000")
	} else{
		library(RColorBrewer)
		colours <- brewer.pal(Ncol, "YlOrRd")
	}

	if(is.null(brks))
	{
		ci <- classInt::classIntervals(x@data[,attr], Ncol, style = "kmeans")
		brks <- ci$brks * 100
		brks <- c(floor(brks[1]), round(brks[-c(1,length(brks))]), ceiling(brks[length(brks)]))/100
	}

	if(class(x) %in% c("SpatialGridDataFrame","SpatialPixelsDataFrame"))
	{
		library(raster)
		x <- raster(x, layer = attr)
		x <- projectRaster(x, crs = "+proj=longlat +datum=NAD83")
	} else {
		stop("x must be SpatialGridDataFrame or SpatialPixelsDataFrame object")
	}
	
	if(is.null(xlim))
	{
		xlim= bbox(x)[1,]
	}
	if(is.null(ylim))
	{
		ylim= bbox(x)[2,]
	}
	
	bgMap <- getTiles(xlim, ylim, zoom, path = "http://a.www.toolserver.org/tiles/osm-no-labels/", maxTiles = 400)
		
	if(!is.null(filename))
	{ 
		if(!length(grep("\\.png$", filename)))
		{
			png(paste(filename, ".png", sep = ""), width = width, height = height)
		}	
	}
	
	plot(1, type = "n", xlim = xlim,ylim=ylim, xlab=NA,ylab=NA,axes=F)

	image(bgMap, add = TRUE)
	
	brkstemp <- brks
	brkstemp[length(brkstemp)] <- brkstemp[length(brkstemp)] + 0.01
	plot(x, breaks = brkstemp, col = paste(colours, substring(hsv(alpha = as.numeric(trans)/100),8), sep = ""), 
		legend = FALSE, horizontal = FALSE, add = TRUE)

	legend("bottomright", fill = colours, cex = label, 
		legend = c(paste("[", brks[1:(length(brks) - 2)], ",", brks[-c(1,length(brks))], ")", sep = ""), 
			paste("[", brks[length(brks) - 1], ",", brks[length(brks)], "]", sep = "")), bg = "white")
	lines(c(-82.21666,-81.97337), c(42.55,42.55), lwd = 4)
	lines(c(-82.21666,-82.21666), c(42.54,42.56), lwd = 2)
	lines(c(-81.97337,-81.97337), c(42.54,42.56), lwd = 2)
	text(-82.095015, 42.52, "20km", cex = label)

	if(!is.null(filename) & devOff)
	{
		dev.off()
	}
}
