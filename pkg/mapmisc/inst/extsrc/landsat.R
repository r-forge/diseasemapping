#+ packages
if(FALSE) {
  install.packages("Pmisc", 
      repos='http://r-forge.r-project.org')
  # also needs R.utils
}
#'

#+ getLandsatTilesFunction
getLandsatTiles = function(x,
    tiles = landsat,
    buffer = c(1,0.6,0.3, 0.6)*50*1000,
    crsOut = 'landsat') {
  
  xExtent = 
      raster(extend(extent(x), buffer), 
          ncol=10, nrow=10, crs=projection(x))
  
  xExtentProj = projectExtent(xExtent, projection(tiles))
  
  tilesCrop = raster::crop(
      tiles, xExtentProj
  )
  
  if(identical(crsOut, 'landsat'))
    crsOut = paste(
      '+proj=lsat +path=',
      tilesCrop$PATH[1],
      ' +lsat=5 +ellps=WGS84 +datum=NAD83', sep='')
  
  tilesProj = spTransform(tiles, crsOut)
  tilesProjCrop = tilesProj[tilesProj$PR %in% 
          tilesCrop$PR,]
  
  tilesProjCrop$pPath = paste('p',  sprintf("%03d", tilesProjCrop$PATH), sep='')
  tilesProjCrop$rRow = paste('r',  sprintf("%03d", tilesProjCrop$ROW), sep='')
  
  tilesProjCrop$pathString = sprintf("%03d", tilesProjCrop$PATH)
  tilesProjCrop$rowString = sprintf("%03d", tilesProjCrop$ROW)
  
  tilesProjCrop$tileString = paste(
      tilesProjCrop$pathString,
      tilesProjCrop$rowString, sep='')
  
  tilesProjCrop$tTile = paste(
      tilesProjCrop$pPath,
      tilesProjCrop$rRow, sep=' ')
  
  tilesProjCrop
}
#'

#+ setup
dataDir = gsub("[[:alnum:]]+$", "landsat", tempdir())
if(!dir.exists(dataDir)) dir.create(dataDir)
#'

#+ downloadLandsatCover
wfiles = Pmisc::downloadIfOld(   'https://landsat.usgs.gov/sites/default/files/documents/wrs2_descending.zip',
    path = dataDir
)

library(raster)
landsatCover = shapefile(grep("shp$", wfiles, value=TRUE))
#'

#+ landsatProj

landsatSub = landsatCover[landsatCover$PATH %in% seq(0,30,by=2),]

landsatCrop = raster::crop(
    landsatSub, extent(-170,170,-80,75)
)
cropImage = mapmisc::openmap(
    landsatCrop
)    
plot(cropImage)
plot(landsatCrop, add=TRUE, border='#FF0000')    
"+proj=ob_tran +o_proj=cc +o_lat_p=50 +o_lon_p=0 +ellps=WGS84 +datum=NAD83"

#  eck2 cc
"+proj=lsat +path=20 +lsat=5 +ellps=WGS84 +datum=NAD83"
"+proj=eck2 +ellps=WGS84 +datum=NAD83" 
landsatTrans =  spTransform(
    landsatCrop, 
    mapmisc::omerc(toLL, angle=15)
)
plot(landsatTrans)

torontoT = spTransform(toLL, landsatTrans@proj4string)

transImage = mapmisc::openmap(
    torontoT, zoom=5, fact=2
)    
mapmisc::map.new(transImage)
plot(transImage,add=TRUE)
plot(landsatTrans, add=TRUE)    

quantile(rgeos::gArea(landsatTrans))/10^12

toTile = landsatTrans[10,]
plot(toTile@polygons[[1]]@Polygons[[1]]@coords)
polygon(toTile@polygons[[1]]@Polygons[[1]]@coords)

#'

#+ landsatToronto
toLL = mapmisc::geocode("Toronto, ON")

toOmerc = spTransform(toLL, 
    mapmisc::omerc(toLL, -17))

landsatTilesO = getLandsatTiles(toOmerc, landsatCover)

plot(landsatTilesO)

landsatCrs = mapmisc::omerc(
    landsatTiles[1,], angle=seq(5,20,by=0.1)
)
landsatTiles = spTransform(landsatTilesO, landsatCrs)
plot(landsatTiles)
#'

#+ torontoMap
torontoTiles = mapmisc::tonerToTrans(
    mapmisc::openmap(landsatTiles, path='stamen-toner')
)

Stiles = seq(1, len=length(landsatTiles))
Scol = mapmisc::col2html(
    RColorBrewer::brewer.pal(
        length(landsatTiles), 
        'Dark2'
    ),
    0.85
)
mapmisc::map.new(landsatTiles, buffer=-2*1000)
plot(torontoTiles, add=TRUE)
plot(landsatTiles, add=TRUE, 
    border=Scol,
    col = Scol, 
    density = 12,
    angle = 180*Stiles/length(Stiles),
    lwd=3
)
mapmisc::legendBreaks("bottomleft", 
    breaks= landsatTiles$tTile, 
    col = Scol)    
#'

#+ selectTiles

years = 2014 
months = c('Feb','Jul')
days = 1:31

forLandsat = list(
    landsat = 'L',
    satellite = 8,
# Sensor: “C” = OLI/TIRS Combined, “E” = ETM+, “T” = TM, “M”= MSS)
    sensor = 'C', 
# GSI = Ground station identifier
    groundStation = '002',
# VV = Archive version number
    version = '00'
)



timeMat = expand.grid(
    tile = landsatTiles$tTile,
    year=years, 
    month=months,
    day = days
)


timeMat$date =strptime(paste(
        timeMat$year, 
        timeMat$month,
        timeMat$day, sep='-'), 
    format='%Y-%b-%d', tz='UTC')
timeMat = timeMat[!is.na(timeMat$date),]


timeMat$julianDay = format(timeMat$date, format='%j')
timeMat$julianTwo = substr(timeMat$julianDay, 1,2)

timeMat = timeMat[!duplicated(
        timeMat[,c('tile','year','julianTwo')]
    ), ]


tileMat = merge(
    timeMat, landsatTiles, 
    by.x='tile', by.y = 'tTile',
    all.x=TRUE, all.y=FALSE
)


tileMat$marker = paste(
    forLandsat$landsat, 
    forLandsat$satellite,
    '/',
    tileMat$pathString,
    '/',
    tileMat$rowString,
    sep=''
)    

tileMat$prefix = paste(
    forLandsat$landsat, 
    forLandsat$satellite,
    '/',
    tileMat$pathString,
    '/',
    tileMat$rowString,
    '/',
    forLandsat$landsat, 
    forLandsat$sensor,
    forLandsat$satellite,
    tileMat$pathString,
    tileMat$rowString,
    tileMat$year, 
    tileMat$julianTwo,
    sep=''
)    


#'

#' https://cran.r-project.org/web/packages/getlandsat/vignettes/getlandsat_vignette.html
#+ getTiles
library("getlandsat")

landsatScenes = mapply(
    getlandsat::lsat_list,
    marker = tileMat$marker,
    prefix = tileMat$prefix,
    MoreArgs = list(max=10000)
)
#'

#+ formatTiles

landsatScenesShort = landsatScenes[
    which(unlist(lapply(landsatScenes, nrow))>0)
]

landsatScenes2 = do.call(
    rbind, landsatScenesShort
)

# ndvi, need band 4 and 5    
landsatScenesForNdvi = landsatScenes2[
    grep("b(4|5)[.]tif$", landsatScenes2$Key, ignore.case=TRUE),
]    

landsatNdvi = data.frame(
    pathString = substr(landsatScenesForNdvi$Key,4,6),
    rowString = substr(landsatScenesForNdvi$Key,8,10),
    year = substr(landsatScenesForNdvi$Key,21,24),
    julianDay = substr(landsatScenesForNdvi$Key,25,27),
    band = gsub("^.*_B|[.]TIF$", "", 
        landsatScenesForNdvi$Key, ignore.case=TRUE),
    file = landsatScenesForNdvi$Key,
    stringsAsFactors=FALSE
)    
landsatNdvi$time= strptime(
    paste(landsatNdvi$year, landsatNdvi$julianDay),
    format = '%Y %j'
)
landsatNdvi$month = format(landsatNdvi$time, '%b')

Simage = sapply(
    landsatNdvi$file,
    getlandsat:::parse_landsat_str
)

landsatNdvi$url = paste(
    "https://s3-us-west-2.amazonaws.com/landsat-pds/L8/",
    landsatNdvi$pathString, '/', 
    landsatNdvi$rowString, '/',
    gsub("_B[[:alnum:]].*[.]tif$", "", 
        basename(landsatNdvi$file), ignore.case=TRUE),
    '/', basename(landsatNdvi$file),
    sep=''
)    

landsatNdvi$localFiles = Pmisc::downloadIfOld(
    landsatNdvi$url, path=dataDir,
    age = '2 years', mode='wb')

#'

#+ oneTile

stuff = stack("C:\\Users\\pbrown\\AppData\\Local\\Temp\\landsat\\LC80170302014048LGN00_B5.TIF")
plot(stuff)


#'
