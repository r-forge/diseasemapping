library(XML)
murder = xmlTreeParse("http://www3.thestar.com/static/googlemaps/homicides.xml" )
canMap = readOGR("/store/patrick/spatialData/gpr_000b11a_e/", "gpr_000b11a_e",
		verbose=F)

murder = murder[[1]][[1]][[3]]

murder = as.data.frame(t(xmlSApply(murder, function(qq) as.numeric(xmlAttrs(qq[["point"]])))))
names(murder) = c("long","lat")
murder= na.omit(murder)
murder = SpatialPoints(murder, proj4string = CRS("+init=epsg:4326"))

murder = murder[!duplicated(murder@coords[,1]),]



torontoCT = readOGR("/store/patrick/spatialData/","gct_000b06a_e")
torontoCT = torontoCT[torontoCT$CMAUID==535,]




#income = read.table("/store/patrick/spatialData/98-316-XWE2011001-401.CSV",#
#		sep=",",header=T,skip=1, quote="\"",as.is=T,nrows=1500000)
#income = income[income$CMA=="Toronto",]

library(diseasemapping)
data(popdata)
toronto = popdata[popdata$CSDNAME=="Toronto",]
murder=murder[!is.na(over(murder, toronto))]

murderUTM = spTransform(murder, 
		CRS("+proj=utm +zone=19 +datum=WGS84 +units=m +no_defs +ellps=WGS84"))

torontoBorder = spTransform(toronto,
		CRS("+proj=utm +zone=19 +datum=WGS84 +units=m +no_defs +ellps=WGS84"))



plot(murderUTM)
plot(spTransform(toronto, murderUTM@proj4string), add=T,border='red')

torTiles = getTiles(toronto@bbox,zoom=11,path="http://tile.openstreetmap.org/")

worldNight = raster("/store/patrick/spatialData/world_avg_dat.tif")
torontoNight = crop(worldNight,extent(toronto@bbox + 0.1*cbind(c(-1,-1),c(1,1))))
down=0.02;right=0.03
torontoNight@extent@xmin = torontoNight@extent@xmin +right
torontoNight@extent@xmax = torontoNight@extent@xmax +right
torontoNight@extent@ymin = torontoNight@extent@ymin -down
torontoNight@extent@ymax = torontoNight@extent@ymax -down


torontoNightUTM = projectRaster(torontoNight, 	
		crs=murderUTM@proj4string)
torontoNightUTM = crop(torontoNightUTM, extent(torontoBorder@bbox))
plot(torontoNightUTM^2)
plot(torontoBorder,add=T)

load("/home/patrick/workspace/Admin/teaching/spatial4/data/tor06.RData")

tor06@data[which(tor06$econ_fam_med_inc==0),"econ_fam_med_inc"] = NA
tor06 = spTransform(tor06, murder@proj4string)

myExtent = extent(tor06@bbox)
therange = apply(tor06@bbox,1,diff)

cellsX = 250
dimX = therange["x"]/cellsX
cellsY = ceiling(cellsX *therange["y"]/therange["x"])
myExtent@ymax = myExtent@ymin + cellsY * dimX

rasterTemplate = raster(myExtent, nrows=cellsY, ncols=cellsX,
		crs= (murder@proj4string))

torontoIncome = rasterize(tor06,rasterTemplate, "econ_fam_med_inc",
		fun=mean, na.rm=T)


tor06$area_sqk=sapply(slot(tor06, "polygons"), slot, "area")/(1000^2)

tor06$pdens_hectare= tor06$pdens = (tor06$m_tot + tor06$f_tot)/
		(sapply(slot(tor06, "polygons"), slot, "area")/10000)

torontoPdens = rasterize(tor06,rasterTemplate, "pdens_hectare",
		fun=mean, na.rm=T)


plot(murderCovariates[["income"]])
plot(murderUTM, add=T)

murder=murderUTM
torontoNight= torontoNightUTM

save(murder,torontoIncome, torontoNight, torontoPdens,
		file="/home/patrick/workspace/diseasemapping/pkg/geostatsp/data/murder.RData",
		compress="xz")

ltLoa = crop(ltLoa, extent(loaloa@bbox))
names(eviLoa) = "evi"


eviLoa = aggregate(eviLoa, 2)
elevationLoa = aggregate(elevationLoa, 2)
save(elevationLoa,eviLoa, loaloa, ltLoa, 
		file="/home/patrick/workspace/diseasemapping/pkg/geostatsp/data/loaloa.RData",
		compress="xz")
