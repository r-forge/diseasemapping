#+ setup

library('mapmisc')

#if(!exists('mapmiscCachePath'))
#  mapmiscCachePath = system.file('extdata', package='mapmisc')

#if(!exists('mapmiscCacheReadOnly'))
#  mapmiscCacheReadOnly = TRUE


options(
#  mapmiscCachePath = mapmiscCachePath,
#  mapmiscCacheReadOnly = mapmiscCacheReadOnly,
  mapmiscVerbose=TRUE)

#getOption("mapmiscCachePath")
#getOption("mapmiscCacheReadOnly")
#'

if(!interactive()) pdf("openmap.pdf")

#+ simplePlot
myraster = rast(ext(8,18,0,10), crs=crsLL)
values(myraster) = seq(0,1,len=ncell(myraster))

myPoints = as.points(myraster)[
  seq(1,ncell(myraster),len=5)]

plot(myraster)
points(myPoints)
#'
#' 
#+ africaPlots 
  
  # utm zone 32
  utmproj = crs("+init=epsg:3064") 
  myrasterUTM = project(myraster, utmproj)
  myPointsUTM = project(myPoints, utmproj)
  plot(myrasterUTM)
  points(myPointsUTM)
  
  myPointsMercator = project(myPoints, crsMerc)
  
  
  myplot = function(first,second=first) {
    par(mar=c(0,0,0,0))
    plot(first)
    plot(mytiles, add=TRUE)
    plot(second,add=TRUE,col='blue')
#	points(mycities,col='red')
#	text(mycities, labels=mycities$name, col='red',pos=4)
    scaleBar(first)
  }
  
  thezoom=6
  
  print(1)

# and if the OpenStreetMap.org web site can be accessed
  
  # raster, result will be in project of the raster (long-lat)
  mytiles = openmap(
    x=extend(myraster,1),
    zoom=thezoom)
#		mycities = GNcities(extend(myraster,1),max=5)
  myplot(myraster, myPoints)

print(2)


  # slash at the end
  mytiles = openmap(extend(myraster,1),zoom=thezoom, 
    path="http://tile.openstreetmap.org/")
#		mycities = GNcities(extend(myraster,1),max=5)
  myplot(myraster, myPoints)
print(3)

  
  # no http at beginning
  mytiles = openmap(extend(myraster,1),path="tile.openstreetmap.org")
#		mycities = GNcities(extend(myraster,1),max=5)
  myplot(myraster, myPoints)
  

print(4)

  
  # extent, tiles will be long-lat
  mytiles = openmap(ext(myraster),zoom=thezoom)
  # cities will be long=lat
#		mycities = GNcities(extent(myraster),max=5,lang="fr")
#		myplot(mycities,myPoints)

print(5)

  
  # give the bbox, long lat
  mytiles = openmap(ext(myraster),zoom=thezoom)
#		mycities = GNcities(bbox(myraster),max=5)
#		myplot(mycities,myPoints)
  

print(6)

  
  # give points, result is CRS of points (long-lat)
  mytiles = openmap(myPoints,zoom=thezoom)
#		mycities = GNcities(myPoints,max=5,lang="es")
  myplot(myPoints)

print(7)

  
  # UTM raster
  mytiles = openmap(myrasterUTM,zoom=thezoom)
#		mycities = GNcities(myrasterUTM,max=5)
  myplot(myrasterUTM, myPointsUTM)

print(8)

  
  # supply a crs
  mytiles = openmap(x=ext(myrasterUTM),zoom=thezoom, 
    crs=crs(myrasterUTM))
#		mycities = GNcities(myrasterUTM,max=5)
  myplot(myrasterUTM, myPointsUTM)

print(9)

  
  # utm points
  mytiles = openmap(myPointsUTM,zoom=thezoom)
#		mycities = GNcities(myPointsUTM,max=5)
  myplot(myPointsUTM)

print(10)

  
  # specify different output crs
  mytiles = openmap(myPointsUTM, crs=crsLL)
#	mycities = GNcities(myPoints,max=5)
  myplot(myPoints)

print(11)

  
  # one point only
  mytiles = openmap(crds(myPoints)[1,], zoom=4, crs=crs(myPoints), buffer=10, fact=2)
  myplot(myPoints)

print(12)  
#'
#'   
#' ams city hall
#+ ams
  cityHall = vect(cbind(4.891111, 52.373056),
    crs=crsLL)
#  cityHall = spTransform(cityHall,CRS("+init=epsg:28992"))
  cityHall = project(cityHall,crs("+init=epsg:32631"))
  mytiles = openmap(cityHall, buffer=50)

  map.new(mytiles)
  plot(mytiles, add=TRUE)
  points(cityHall, pch=3, col='blue',cex=4)
  scaleBar(mytiles, 'topleft', bty='n', col='red')
#'  
  
  if(!interactive()) dev.off()
