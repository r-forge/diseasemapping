#+ setup
dataDir = '/store/patrick/spatialData'
library("rgdal")
library('raster')
#install.packages("MODIS", repos="http://R-Forge.R-project.org")
library("MODIS")
#'

#+ map of the world
worldUrl = "http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/cultural/ne_10m_admin_0_countries_lakes.zip"
worldZip = file.path(dataDir,"foo.zip")
if(!file.exists(worldZip))
	download.file(worldUrl, worldZip)
unzip(worldZip,
		exdir=dataDir)
unzip(worldZip,list=TRUE)
worldShapeFile=grep("\\.shp$", unzip(worldZip,list=TRUE)$Name,value=TRUE)
worldShapeFile = gsub("\\.shp$","",worldShapeFile)
world = readOGR(dataDir,worldShapeFile,
		stringsAsFactors=FALSE,verbose=FALSE)
indiaChina = world[world$NAME%in%c("India","China"),]
thebox = extent(indiaChina)
#'

#' sulpher dioxide
#' http://mirador.gsfc.nasa.gov/cgi-bin/mirador/presentNavigation.pl?tree=project&dataset=OMSO2e.003&project=OMI&dataGroup=L3_V003&version=003
#+ downloadSO2Files 
library(RCurl)
bob=getURL('ftp://acdisc.gsfc.nasa.gov/data/s4pa///Aura_OMI_Level3/OMSO2e.003/2010/',
		ftp.use.epsv=TRUE,
		dirlistonly = TRUE)
bob=strsplit(bob, '\n')[[1]] 
feb = sort(grep("^[[:print:]]+2010m02[[:print:]]+he5$",bob,value=T))[1:5]


for(D in feb) {
	thisfile =file.path(dataDir,D) 
	if(!file.exists(thisfile)) {
		download.file(
				paste(
						'ftp://acdisc.gsfc.nasa.gov/data/s4pa///Aura_OMI_Level3/OMSO2e.003/2010/',
						D,sep=""),
			thisfile)
	}
}
theFiles = file.path(dataDir, feb)
#'



#+ readSO2files,cache=TRUE
theLayers = MODIS::getSds(theFiles[12])

x = stack(theLayers$SDS4gdal)
extent(x) = extent(-180,180,-90,90)
projection(x) = CRS("+init=epsg:4326")
x = crop(x, thebox)

x2 =  reclassify(x,
		matrix(c(-Inf, -100, NA),ncol=3,byrow=T))
names(x2) = names(x)

Dvar = "ColumnAmountSO2_PBL"
x2[[Dvar]] =   reclassify(x2[[Dvar]] ,
		matrix(c(10, Inf, NA),ncol=3,byrow=T))

#'

#+ downloadMap
library('mapmisc')
bgmap = openmap(x2)
#'

#+ aMap
thevar = c('ColumnAmountO3','ColumnAmountSO2_PBL')
par(mfrow=c(1,2))
for(Dvar in thevar) {
theColours = colourScale(x2[[Dvar]], style='quantile',breaks=7, dec=1,opacity=0.5)

map.new(x2)
plot(bgmap,add=TRUE)
plot(x2[[Dvar]], add=TRUE, col=theColours$colOpacity, breaks=theColours$breaks, legend=FALSE)
legendBreaks("bottomright", theColours)
}
#'



http://disc.sci.gsfc.nasa.gov/measures/documentation/seawifs-deepblue-level-3-data-fields

#+ moreStuff
thisfile =file.path(dataDir,"bc.h5") 
if(!file.exists(thisfile)) {
	download.file("ftp://measures.gsfc.nasa.gov/data/s4pa//DeepBlueSeaWiFS_Level3/SWDB_L305.004/2010/DeepBlue-SeaWiFS-0.5_L3_20101128_v004-20130604T135711Z.h5",
			thisfile)
}
theLayers = MODIS::getSds(thisfile)
x = stack(theLayers$SDS4gdal)
extent(x) = extent(-180,180,-90,90)
projection(x) = CRS("+init=epsg:4326")
x = crop(x, thebox)

x2 =  reclassify(x,
		matrix(c(-Inf, -100, NA),ncol=3,byrow=T))
names(x2) = names(x)

D='angstrom_exponent_ocean'
plot(x2[[D]])
plot(world,add=TRUE)


#'