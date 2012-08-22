plotWithTiles = function(x, attr = 1, breaks = NULL, prob = FALSE, 
	Ncol = 5, trans = 50, zoom = 4, xlim = NULL, ylim = NULL,
	filename = NULL, devOff = TRUE, width = 6, height = 6, label = 1,
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
		
		D=1
		breaks = signif(ci$brks, D)
		
		
		while(length(unique(breaks)) < length(breaks)) {
			D = D+1
			breaks = signif(ci$brks, D)
		}
			
		themin = min(x@data[,attr], na.rm=T)
		themax = max(x@data[,attr], na.rm=T)
		
		if(themin < breaks[1]) {
			 breaks[1] = breaks[1] - 10^(-(D+1))
		}
		
		if(themax > breaks[length(breaks)]) {
			breaks[length(breaks)] = breaks[length(breaks)] + 10^(-D+1)
		}
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
		if(class(x)=='raster'|class(x)=='RasterLayer') {
		xlim= apply(as.matrix(x), 2, function(qq) all(is.na(qq)))
		xlim= seq(bbox(x)['s1','min'], bbox(x)['s1', 'max'], len=x@ncols)[!xlim]
	    xlim = range(xlim)
		} else
			xlim = bbox(x)[1,] 
	}
	if(is.null(ylim))
	{
		if(class(x)=='raster'|class(x)=='RasterLayer') {
			ylim= apply(as.matrix(x), 1, function(qq) all(is.na(qq)))
		ylim= seq(bbox(x)['s2','min'], bbox(x)['s2', 'max'], len=x@nrows)[!ylim]
	    ylim = range(ylim)
	} else
		ylim = bbox(x)[2,] 
	
	}
	
	bgMap <- getTiles(xlim, ylim, zoom, path = tilepath, maxTiles = 400)
		
	if(!is.null(filename))
	{ 
		if(!length(grep("\\.png$", filename)))
		{
			filename = paste(filename, ".png", sep = "")
		}	
		png(filename, width = width, height = height,units='in',res=300)

	}

	dummy=SpatialPoints(cbind(xlim, ylim), proj4string=CRS("+proj=longlat +datum=NAD83"))
	par(mar=c(0,0,0,0))
	plot(dummy,  xlim = xlim,ylim=ylim, xlab=NA,ylab=NA,axes=F)

	image(bgMap, add = TRUE)

	if(!is.null(attr)) {
	plot(x, breaks = breaks, col = paste(colours, substring(hsv(alpha = as.numeric(trans)/100),8), sep = ""), 
		legend = FALSE, horizontal = FALSE, add = TRUE)

	
	
legend('topleft', pch=rep(NA, length(breaks)), legend=breaks, bg='white',
		adj=c(0,-0.4), cex=label)

legend("topleft", pch=15, col =  colours, pt.cex = 2.75*label, cex=label, 
			legend=rep('    ', length(colours)), bty='n')
	
	
	} else {
		plot(x, add=T, col=colours)		
	}

	if(!is.null(filename) & devOff)
	{
		dev.off()
	}
}
