plotWithTiles = function(x, attr=1, brks = NULL, prob = FALSE, 
	Ncol = 5, trans = 50, zoom = 4, xlim=NULL,ylim=NULL,
	filename = NULL, devOff = TRUE, width = 1200, height = 1200
)
{
	library(classInt)
	library(webmaps)
	library(spdep)
	
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
		ci <- classIntervals(x@data[,attr], Ncol, style = "kmeans")
		brks <- ci$brks * 100
		brks <- c(floor(brks[1]), round(brks[-c(1,length(brks))]), ceiling(brks[length(brks)]))/100
	} else{
		ci <- classIntervals(x@data[,attr], length(brks), fixedBreaks = brks, 
				style = "fixed")
	}


	if(class(x)%in%c("SpatialGridDataFrame","SpatialPixelsDataFrame")) {

			library(rgdal)
		#x@data = x@data[,attr,drop=F] 

		mytempfile = paste(tempfile(), ".tif",sep="")
		mytempfile2 = paste(tempfile(), ".tif",sep="")
		
		writeGDAL(x[attr], mytempfile)
	
		
		thebbox=x@bbox
		
		
		thebbox = 	SpatialPoints(
				expand.grid(thebbox[1,],thebbox[2,]),
				 proj4string=x@proj4string)
		
		 thebbox = spTransform(thebbox, CRS("+proj=longlat +datum=NAD83"))
		 
		 newX = sort(thebbox@coords[,1])
		 newX = newX[-c(1,length(newX))]
		 newY = sort(thebbox@coords[,2])
		 newY = newY[-c(1,length(newY))]
		 
		 if(.Platform$OS.type =="unix") {
		     
		system(
				paste("gdalwarp -s_srs '", x@proj4string@projargs, 
						"' -t_srs '+proj=longlat +datum=NAD83' -te",
						newX[1], newY[1], newX[2], newY[2],
    					mytempfile, mytempfile2, sep = " ")
		)
		
		
		} else {
		# put windows code in here.
			warning("can't do this on windows yet")
		}

		x <- readGDAL(mytempfile2)
		file.remove(mytempfile)
		file.remove(mytempfile2)
		
	} else {
		x <- spTransform(x, CRS("+proj=longlat +datum=NAD83"))
	}
	

	if(is.null(xlim))
		xlim= x@bbox[1,] 
	if(is.null(ylim))
		ylim= x@bbox[2,]
	
	bgMap <- getTiles(xlim, ylim, zoom, path = "http://tile.openstreetmap.org/", 
				maxTiles = 200)
	
	
	if(!is.null(filename))
	{ 
		png(paste(filename, ".png", sep = ""), width = width, height = height)
	}
	
	
	
	
	plot(x, xlim = xlim,ylim=ylim, type="n",xlab="longitude",ylab="latitude",
			border="white",col="white")

	image(bgMap, add = TRUE)
	
	if(class(x) == "SpatialGridDataFrame")
	{
		image(x, "band1", 
				col=paste(colours, 
						substring(hsv(alpha = as.numeric(trans)/100),8), sep = ""), 
				breaks = brks, add = TRUE)
	} else {
			colFac <- findColours(ci, colours)
			colFac[!is.na(colFac)] = paste(colFac[!is.na(colFac)], trans,sep="")
			
			plot(x, col = colFac, add=T)
	}
	legend("bottomright", fill = colours, legend = brks[-1], bg = "white")

	if(!is.null(filename) & devOff )
	{
		dev.off()
	}
}


testPlotSurface= function() {
	
	library(diseasemapping)
	data(popdata)
	proj4string(popdata) = CRS("+proj=longlat +datum=NAD83")
	install.packages("pixmap")
	install.packages("webmaps", repos="http://R-Forge.R-project.org")
	plotWithTiles(popdata, 0, zoom=6,
    		trans=40,
	xlim=c(-81.95318, -78.45677),
	ylim=c(42.05731, 44.94908)	
)


yseq = seq(5500000, len=500,by=1000)
xseq= seq(0,len=500,by=1000)
			
stuff = SpatialPixelsDataFrame(expand.grid(xseq,yseq),
		proj4string=CRS("+init=epsg:26917"), 
		data=data.frame(z=rnorm(length(xseq)*length(yseq))))


plotWithTiles(stuff,zoom=6)

}
