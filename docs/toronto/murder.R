library(XML)
murder = xmlTreeParse("http://www3.thestar.com/static/googlemaps/homicides.xml" )
canMap = readOGR("/store/patrick/spatialData/gpr_000b11a_e/", "gpr_000b11a_e",
		verbose=F)

murder = murder[[1]][[1]][[3]]

murder = as.data.frame(t(xmlSApply(murder, function(qq) as.numeric(xmlAttrs(qq[["point"]])))))
names(murder) = c("long","lat")
murder= na.omit(murder)
murder = SpatialPoints(murder, proj4string = CRS("+init=epsg:4326"))




torontoCT = readOGR("/store/patrick/spatialData/","gct_000b06a_e")
torontoCT = torontoCT[torontoCT$CMAUID==535,]

#income = read.table("/store/patrick/spatialData/98-316-XWE2011001-401.CSV",#
#		sep=",",header=T,skip=1, quote="\"",as.is=T,nrows=1500000)
#income = income[income$CMA=="Toronto",]

library(diseasemapping)
toronto = popdata[popdata$CSDNAME=="Toronto",]
murder=murder[!is.na(overlay(murder, toronto))]

murderUTM = spTransform(murder, 
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
plot(torontoNightUTM)

load("../data/tor06.RData")
tor06@data[which(tor06$econ_fam_med_inc==0),"econ_fam_med_inc"] = NA

murderCovariates = stackRasterList(list(night=torontoNightUTM, income=tor06[,"econ_fam_med_inc"]),
		template=torontoNightUTM)

plot(murderCovariates[["income"]])
plot(murderUTM, add=T)

murder=murderUTM
save(murder,murderCovariates, 
		file="/home/patrick/workspace/diseasemapping/pkg/geostatsp/data/murder.RData",
		compress="xz")

