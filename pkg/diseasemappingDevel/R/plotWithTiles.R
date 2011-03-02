plotWithTiles = function(x, attr = 1, brks = NULL, prob = FALSE, 
	Ncol = 5, trans = 50, zoom = 4, xlim = NULL, ylim = NULL,
	filename = NULL, devOff = TRUE, width = 600, height = 600
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

	if(!is.na(attr)) 
	{
		if(is.null(brks))
		{
			ci <- classInt::classIntervals(x@data[,attr], Ncol, style = "kmeans")
			brks <- ci$brks * 100
			brks <- c(floor(brks[1]), round(brks[-c(1,length(brks))]), ceiling(brks[length(brks)]))/100
		} else {
			breaksQQQ <<- brks
			ci <- classInt::classIntervals(x@data[,attr], n = length(brks) - 1, fixedBreaks = breaksQQQ, style = "fixed")
			rm(breaksQQQ, pos = 1)
		}
	}

	if(class(x)%in%c("SpatialGridDataFrame","SpatialPixelsDataFrame")) 
	{
		library(rgdal)
		#x@data = x@data[,attr,drop=F] 


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

		 	mytempfile = paste(tempfile(), ".tif",sep="")
			mytempfile2 = paste(tempfile(), ".tif",sep="")
			writeGDAL(x[attr], mytempfile)
		 
		system(
				paste("gdalwarp -s_srs '", x@proj4string@projargs, 
						"' -t_srs '+proj=longlat +datum=NAD83' -te",
						newX[1], newY[1], newX[2], newY[2],
    					mytempfile, mytempfile2, sep = " ")
		)
		
		
		} else {
	
			mytempfile = "temporaryfile1.tif"
			mytempfile2 = "temporaryfile2.tif"
			writeGDAL(x[attr], mytempfile)
			
			system(
			paste('cmd /c c:\\OSGeo4W\\bin\\gdalwarp.exe', 
				mytempfile, mytempfile2, 
				'-s_srs "', x@proj4string@projargs, 
				'" -t_srs "+proj=longlat +datum=NAD83" -te',
				newX[1], newY[1], newX[2], newY[2], 
				sep = " "))
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
	
	bgMap <- getTiles(xlim, ylim, zoom, path = "http://tile.openstreetmap.org/", maxTiles = 200)
	
	
	if(!is.null(filename))
	{ 
		if(!length(grep("\\.png$", filename)))
			png(paste(filename, ".png", sep = ""), width = width, height = height)
	}
	
	plot(x, lty = 0, xlim = xlim,ylim=ylim, xlab="longitude",ylab="latitude",col=NA,axes=F)

	image(bgMap, add = TRUE)
	
	if(class(x) == "SpatialGridDataFrame")
	{
		image(x, "band1", 
				col=paste(colours, 
						substring(hsv(alpha = as.numeric(trans)/100),8), sep = ""), 
				breaks = brks, add = TRUE)
	} else {
		if(!is.na(attr)) 
		{
			colFac <- findColours(ci, colours)
			colFac[!is.na(colFac)] = paste(colFac[!is.na(colFac)], trans,sep="")
			plot(x, lty = 0, col = colFac, add = TRUE)
		} else {
			colFac = "black"
			plot(x, lty = 1, col = colFac, add = TRUE)
		}  

		
	}
	
	if(!is.na(attr)) 
	{
		legend("bottomright", fill = colours, legend = brks[-1], bg = "white")
    }

	if(!is.null(filename) & devOff)
	{
		dev.off()
	}
}
