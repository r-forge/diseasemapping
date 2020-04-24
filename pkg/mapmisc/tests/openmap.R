#+ setup
library('mapmisc')

if(!exists('mapmiscCachePath'))
  mapmiscCachePath = system.file('extdata', package='mapmisc')

if(!exists('mapmiscCacheReadOnly'))
  mapmiscCacheReadOnly = TRUE


options(
  mapmiscCachePath = mapmiscCachePath,
  mapmiscCacheReadOnly = mapmiscCacheReadOnly,
  mapmiscVerbose=TRUE)

getOption("mapmiscCachePath")
getOption("mapmiscCacheReadOnly")
#'


#+ simplePlot
myraster = raster(matrix(0,10,10),xmn=8,xmx=18,ymn=0,ymx=10, 
  crs=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
)
values(myraster) = seq(0,1,len=ncell(myraster))

myPoints = SpatialPoints(myraster, 
  proj4string=CRS(proj4string(myraster)))[
  seq(1,ncell(myraster),len=5)]

plot(myraster)
points(myPoints)
#'
#' 
#' only do the following if rgdal is available
#+ africaPlots 
if(require('rgdal', quietly=TRUE)) {
  
  # utm zone 32
  utmproj = CRS("+init=epsg:3064") 
  myrasterUTM = projectRaster(myraster, crs=utmproj)
  myPointsUTM = spTransform(myPoints, utmproj)
  plot(myrasterUTM)
  points(myPointsUTM)
  
  myPointsMercator = spTransform(myPoints, 
    crsMerc)
  
  
  myplot = function(first,second=first) {
    par(mar=c(0,0,0,0))
    plot(first)
    plot(mytiles, add=TRUE)
    plot(second,add=TRUE,col='blue')
#	points(mycities,col='red')
#	text(mycities, labels=mycities$name, col='red',pos=4)
    scaleBar(first)
  }
}  
  
  thezoom=6
  
# only do the following if running unix (because nsl is available)
# and if the OpenStreetMap.org web site can be accessed
if(require('rgdal', quietly=TRUE)) {
  
  # raster, result will be in project of the raster (long-lat)
  mytiles = openmap(
    x=extend(myraster,1),
    zoom=thezoom)
#		mycities = GNcities(extend(myraster,1),max=5)
  myplot(myraster, myPoints)
}  
if(require('rgdal', quietly=TRUE)) {

  # slash at the end
  mytiles = openmap(extend(myraster,1),zoom=thezoom, 
    path="http://tile.openstreetmap.org/")
#		mycities = GNcities(extend(myraster,1),max=5)
  myplot(myraster, myPoints)
}  
if(require('rgdal', quietly=TRUE)) {
  
  # no http at beginning
  mytiles = openmap(extend(myraster,1),path="tile.openstreetmap.org")
#		mycities = GNcities(extend(myraster,1),max=5)
  myplot(myraster, myPoints)
  
}  
if(require('rgdal', quietly=TRUE)) {
  
  # extent, tiles will be long-lat
  mytiles = openmap(extent(myraster),zoom=thezoom)
  # cities will be long=lat
#		mycities = GNcities(extent(myraster),max=5,lang="fr")
#		myplot(mycities,myPoints)
}  
if(require('rgdal', quietly=TRUE)) {
  
  # give the bbox, long lat
  mytiles = openmap(bbox(myraster),zoom=thezoom)
#		mycities = GNcities(bbox(myraster),max=5)
#		myplot(mycities,myPoints)
  
}  
if(require('rgdal', quietly=TRUE)) {
  
  # give points, result is CRS of points (long-lat)
  mytiles = openmap(myPoints,zoom=thezoom)
#		mycities = GNcities(myPoints,max=5,lang="es")
  myplot(myPoints)
}  
if(require('rgdal', quietly=TRUE)) {
  
  # UTM raster
  mytiles = openmap(myrasterUTM,zoom=thezoom)
#		mycities = GNcities(myrasterUTM,max=5)
  myplot(myrasterUTM, myPointsUTM)
}  
if(require('rgdal', quietly=TRUE)) {
  
  # supply a crs
  mytiles = openmap(x=extent(myrasterUTM),zoom=thezoom, 
    crs=proj4string(myrasterUTM))
#		mycities = GNcities(myrasterUTM,max=5)
  myplot(myrasterUTM, myPointsUTM)
}  
if(require('rgdal', quietly=TRUE)) {
  
  # utm points
  mytiles = openmap(myPointsUTM,zoom=thezoom)
#		mycities = GNcities(myPointsUTM,max=5)
  myplot(myPointsUTM)
}  
if(require('rgdal', quietly=TRUE)) {
  
  # specify different output crs
  mytiles = openmap(myPointsUTM, crs=CRS("+init=epsg:4326"))
#	mycities = GNcities(myPoints,max=5)
  myplot(myPoints)
}  
if(require('rgdal', quietly=TRUE)) {
  
  # one point only
  mytiles = openmap(coordinates(myPoints)[1,], zoom=4)
  myplot(myPoints)
  } # end have rgdal
#'
#'   
#' ams city hall
#+ ams
  cityHall = SpatialPoints(cbind( 4.891111, 52.373056), proj4string=crsLL)
#  cityHall = spTransform(cityHall,CRS("+init=epsg:28992"))
if(require('rgdal', quietly=TRUE)) {
  cityHall = spTransform(cityHall,CRS("+init=epsg:32631"))
  mytiles = openmap(cityHall, buffer=50)

  map.new(mytiles)
  plot(mytiles, add=TRUE)
  points(cityHall, pch=3, col='blue',cex=4)
  scaleBar(mytiles, 'topleft')
} # end have rgdal
#'  
  

