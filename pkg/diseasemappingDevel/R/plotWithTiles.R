plotWithTiles = function(x, attr = 1, breaks = NULL, prob = FALSE, 
	Ncol = 5, trans = 50, zoom = 4, xlim = NULL, ylim = NULL,
	filename = NULL, devOff = TRUE, width = 600, height = 600, label = 1,
	withNames = T 
)
{
	library(classInt)
	library(webmaps)

	
	if(withNames){
		tilepath = "http://tile.openstreetmap.org/"
	} else {
		tilepath = "http://a.www.toolserver.org/tiles/osm-no-labels/"
	}
	
	
	if(is.null(breaks) & !is.null(attr))
	{		
		ci <- classInt::classIntervals(x@data[,attr], Ncol, style = "kmeans")
		
		D=2
		breaks = signif(ci$brks, D)
		
		
		while(length(unique(breaks)) < length(breaks)) {
			D = D+1
			breaks = signif(ci$brks, D)
		}
			
		theunder = min(x@data[,attr], na.rm=T) - breaks[1]
		theover = max(x@data[,attr], na.rm=T) - breaks[length(breaks)]
 
		if(theunder > 0)
			breaks[1] = signif(ci$brks[1] - theunder, D)
		if(theover > 0)
			breaks[length(breaks)] = signif(ci$brks[length(breaks)] + theover, D)
 
	}
	if(is.character(Ncol)){
		colours=Ncol
	} else {
		library(RColorBrewer)
		colours <- rev(brewer.pal(Ncol, "RdYlGn"))
	}
	
	if(prob == TRUE)
	{
		breaks <- c(0,0.2,0.5,0.8,0.95,1)
		colours <- c("#00FF00","#CCEE00","#FFEE44", "#FF6600","#FF0000")
	} 
	
	if(class(x) %in% c("SpatialGridDataFrame","SpatialPixelsDataFrame"))
	{
		library(raster)
		x <- raster(x, layer = attr)

	} 

	if (! class(x) %in% c("RasterLayer","raster", "SpatialPointsDataFrame")  ){
		stop("x must be raster or SpatialGridDataFrame or SpatialPixelsDataFrame or SpatialPointsDataFrame object")
	}

	if(class(x) %in% c("raster","RasterLayer") ){
		x <- projectRaster(x, crs = "+proj=longlat +datum=NAD83")
	} else {
		x = spTransform(x, CRS("+proj=longlat +datum=NAD83"))
	}
	
	if(is.null(xlim))
	{
		xlim= bbox(x)[1,]
	}
	if(is.null(ylim))
	{
		ylim= bbox(x)[2,]
	}
	
	bgMap <- getTiles(xlim, ylim, zoom, path = tilepath, maxTiles = 400)
		
	if(!is.null(filename))
	{ 
		if(!length(grep("\\.png$", filename)))
		{
			png(paste(filename, ".png", sep = ""), width = width, height = height)
		}	
	}

	dummy=SpatialPoints(cbind(xlim, ylim), proj4string=CRS("+proj=longlat +datum=NAD83"))
	plot(dummy,  xlim = xlim,ylim=ylim, xlab=NA,ylab=NA,axes=F)

	image(bgMap, add = TRUE)

	if(!is.null(attr)) {
	plot(x, breaks = breaks, col = paste(colours, substring(hsv(alpha = as.numeric(trans)/100),8), sep = ""), 
		legend = FALSE, horizontal = FALSE, add = TRUE)

	legend("topleft", fill = colours, cex = label, 
		legend = c(paste("[", breaks[1:(length(breaks) - 2)], ",", breaks[-c(1,length(breaks))], ")", sep = ""), 
			paste("[", breaks[length(breaks) - 1], ",", breaks[length(breaks)], "]", sep = "")), bg = "white")
	} else {
		plot(x, add=T, col=colours)		
	}

	if(!is.null(filename) & devOff)
	{
		dev.off()
	}
}
