
#+ setup

source(system.file(file.path('extsrc', 'landsat.R'),
	package = 'mapmisc'))

dataDir = '/store/patrick/landsat'

if(!dir.exists(dataDir)) dir.create(dataDir)
#'

#+ downloadLandsatCover
library(raster)
landsatCover = getLandsatCoverage(dataDir)
landsatScenes = getLandsatSceneList(dataDir)
#'


#+ landsatToronto
toLL = SpatialPoints(
  cbind(-79.3916043, 43.7069564),
  proj4string=mapmisc::crsLL
)

landsatTiles = getLandsatTiles(
  x=toLL, 
  tiles = landsatCover, 
  buffer=0.25)

torontoTiles = mapmisc::tonerToTrans(
  mapmisc::openmap(landsatTiles, path='stamen-toner')
)
#'

#+ torontoLandastPlot
Stiles = seq(1, len=length(landsatTiles))
Scol = mapmisc::col2html(
  RColorBrewer::brewer.pal(
    length(landsatTiles), 
    'Dark2'
  ),
  0.85
)
mapmisc::map.new(landsatTiles)
plot(torontoTiles, add=TRUE)
plot(landsatTiles, add=TRUE, 
  border=Scol,
  col = Scol, 
  density = 12,
  angle = 180*Stiles/length(Stiles),
  lwd=3
)
plot(extent(attributes(landsatTiles)$landsat), 
  add=TRUE, col='grey', lwd=5, lty=3)
points(spTransform(toLL, projection(landsatTiles)),
  col='grey', cex=3, pch=16)
mapmisc::legendBreaks("bottomleft", 
  breaks= landsatTiles$tTile, 
  col = Scol)    
#'

#+ selectTiles
scenesWinter = getLandsatScenes(
  timeRange = ISOdate(
    2013, c(10,12), c(15,15), tz='utc'
  ),
  tiles = landsatTiles,
  sceneList = landsatScenes
)
#'

#+ downloadTiles
landsatWinter = downloadLandsatTiles(
  scenes=scenesWinter, band = c(4,5,'QA'),
  path = dataDir
)  
#'

#+ crop
toOmerc = spTransform(toLL, omerc(toLL, angle=-17))

torontoRegion = raster(
  extend(
    extent(toOmerc),
    c(5000,3000,5000,1000)
  ), 
  res = res(landsatWinter[[1]]),
  crs=projection(toOmerc)
)


landsatCropFiles = gsub('_B[[:digit:]].TIF', '_crop.grd', 
  unlist(lapply(landsatWinter, function(x) filename(x[[1]]))))

for(Dfile in names(landsatWinter)) {
  stuff = projectRaster(
    from = landsatWinter[[Dfile]],
    to = torontoRegion,
    method = 'ngb',
    filename = landsatCropFiles[Dfile],
    overwrite = file.exists(landsatCropFiles[Dfile])
  )    
}

landsatWinterStack = stack(landsatCropFiles)  
#'



#+ moreMapTiles
torontoTiles2 = mapmisc::tonerToTrans(
  mapmisc::openmap(
    landsatWinterStack, 
    path='stamen-toner', 
    fact = 1.5, maxTiles=10)
)

#'

#+ cropMap

febCol = mapmisc::colourScale(
  landsatWinterStack[[1]], breaks=12, 
  dec=-log10(500), col='Spectral',
  style= 'equal', rev=TRUE, transform=-1
)


mapmisc::map.new(torontoRegion)
plot(landsatWinterStack[[1]], col=febCol$col, breaks=febCol$breaks, 
  legend=FALSE, add=TRUE)
plot(torontoTiles2, add=TRUE)
mapmisc::legendBreaks("topleft", febCol)

#'

#+ calcNdvi
ndviBrick = landsatNdvi(landsatWinterStack, path=dataDir)
#'




#+ plotNdviWinter
myCol = mapmisc::colourScale(
  ndviBrick[['mean']], 
  breaks=12, 
  col='BrBG',
  style= 'equal',
  dec=-log10(0.05)
)

mapmisc::map.new(ndviBrick)
plot(ndviBrick[['mean']], col=myCol$col, 
  breaks=myCol$breaks, legend=FALSE, add=TRUE)
plot(torontoTiles2, add=TRUE)
mapmisc::legendBreaks("topleft", myCol)

plot(ndviBrick[['N']])
plot(torontoTiles2, add=TRUE)
#'
