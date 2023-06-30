
getRowCol <- function(extMerc,
  zoom, 
  rasterSphere = .getRasterMerc(zoom)){
  
  Sy=terra::rowFromY(rasterSphere, c(
      terra::ymax(extMerc), terra::ymin(extMerc)
    ))
  Sx=terra::colFromX(rasterSphere, c(
      terra::xmin(extMerc), terra::xmax(extMerc)
    ))
  
  
  list(
    col = seq(Sx[1],Sx[2]),
    row = seq(Sy[1],Sy[2])
  ) 
}


getRasterNrcanDontNeed = function(zoom) {
  
  # raster for maps from 
  # http://geoappext.nrcan.gc.ca/arcgis/rest/services/BaseMaps/CBMT_CBCT_GEOM_3978/MapServer
  
  origin = c(-3.46558E7, 3.931E7) 
  nrCrs = terra::crs(
  "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=49 +lon_0=-95 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs") 
  
  
  nrRes = c(38364.660062653464,
    22489.62831258996,
    13229.193125052918,
    7937.5158750317505,
    4630.2175937685215,
    2645.8386250105837,
    1587.5031750063501,
    926.0435187537042,
    529.1677250021168,
    317.50063500127004,
    185.20870375074085,
    111.12522225044451,
    66.1459656252646,
    38.36466006265346,
    22.48962831258996,
    13.229193125052918,
    7.9375158750317505,
    4.6302175937685215
  )
  
  eps=0.1
  worldLL = terra::rast(
    terra::ext(-180+eps,180-eps,-90+eps,90-eps), 
    crs=crsLL,
    res=0.01
  )
  
  worldNrcan = terra::project(worldLL, nrCrs, res=1000)
  
  nrcanExtent = terra::ext(
    -4282638.06150141,
    4852210.1755664,
    -5153821.09213678,
    4659267.000000001)
  
  nrcanExtent = terra::ext(
    -7786476.885838887,
    7148753.233541353,
    -927807.6591668031,
    7928343.534071138
  )
  
  nrcanRaster = terra::rast(
    nrcanExtent,
    nrow = 2^(1+zoom),
    ncol = 2^(1+zoom),
    crs = nrCrs
  )		
  
  N = 2^zoom 
  
  terra::rast(extentMerc, nrows = N, ncols=N, crs=crsMerc)
}

getRowColNrcanDontNeed <- function(
  extMerc,
  zoom, 
  rasterSphere = .getRasterMerc(zoom)){
  
  
  Sy=terra::rowFromY(rasterSphere, c(
      terra::ymax(extMerc), terra::ymin(extMerc)
    ))
  Sx=terra::colFromX(rasterSphere, c(
      terra::xmin(extMerc), terra::xmax(extMerc)
    ))
  
  
  list(
    col = seq(Sx[1],Sx[2]),
    row = seq(Sy[1],Sy[2])
  ) 
}


nTilesMerc <- function(extMerc,zoom){
  
  SrowCol = terra::getRowCol(extMerc, zoom=zoom)
  
  length(SrowCol[[1]])*length(SrowCol[[2]])
  
}

if(F){
 x = vect(as.matrix(expand.grid(seq(-5e6,-1e6,len=100), seq(-3e6,3e6,len=100))), crs='EPSG:3573')
 xLL = project(x, crsLL)
 zoom=4
 theExt = ext(-4e6,-1e6,-3e6,0e6)
 xx = getTilesFromPoints(
  outraster=rast(theExt, res = (xmax(theExt)-xmin(theExt))/460, crs=crs('EPSG:3573')),
  path="https://stamen-tiles-d.a.ssl.fastly.net/toner/",  
#  path="https://tiles.stadiamaps.com/styles/stamen_watercolor/", suffix='.jpg',
  zoom=3, verbose=TRUE)

}


getTilesFromPoints = function(
  outraster, 
  zoom=1, 
  path="http://tile.openstreetmap.org/",
  cachePath='.',
  cacheDir= make.names(gsub(
          "^http.*//([[:alpha:]][.])*((tile|basemap)s?[.][[:digit:]]?)?(openstreetmap[.])?|[[:punct:]]$", 
          "", path)),
  verbose=FALSE, suffix = '.png',
  tileNames = 'zxy'){
  
  NrowsPerCycle = 20 # number of rows of out raster to process simultaneously

  cacheDir2 = file.path(cachePath, cacheDir)
  
  rasterSphere = .getRasterMerc(zoom)  
  
  samplePoints = rast(ext(outraster), res= (xmax(x)-xmin(x))/100, crs=crs(outraster))
  samplePoints = vect(xyFromCell(samplePoints, 1:terra::ncell(samplePoints)), crs=crs(outraster))
  xMerc = terra::project(samplePoints, crsMerc)

  values(rasterSphere) = NA
  SrowColFull = terra::extract(rasterSphere, xMerc, cells=TRUE)[,'cell']
  SrowColFull = cbind(cell = SrowColFull, row = rowFromCell(rasterSphere, SrowColFull), col = colFromCell(rasterSphere, SrowColFull))

  SrowColFull = SrowColFull[!is.na(SrowColFull[,1]), ]
  SrowCol = SrowColFull[!duplicated(SrowColFull[,c('row','col')]), ]#c('row','col')]
  SrowCol = SrowCol[order(SrowCol[,'col'], SrowCol[,'row']),]
  SrowCol = as.data.frame(SrowCol)
#  SrowCol$newBlock = !c(FALSE, diff(SrowCol[,'col'])==0 & diff(SrowCol[,'row'])==1)
#  SrowCol$block = cumsum(SrowCol$newBlock)


      if(tileNames == 'zxy') {
        SrowCol$cache = file.path(cacheDir2, zoom, SrowCol[,'col']-1)
        SrowCol$path = paste(path,zoom,'/',SrowCol[,'col']-1,'/',sep='')
        SrowCol$tile = paste(SrowCol[,'row']-1, suffix,sep='')
        SrowCol$url = paste(SrowCol$path, SrowCol$tile, sep='')
      } else if (tileNames == 'zxy=') {
        SrowCol$cache = cacheDir2
        SrowCol$tile = paste('z=',zoom, '&x=', SrowCol[,'col']-1, '&y=',SrowCol[,'row']-1,suffix, sep='')
        SrowCol$url = paste(path, SrowCol$tile, sep='')       
      } else if (tileNames == 'xyz=') {
        SrowCol$cache = cacheDir2
        SrowCol$tile = paste('x=',SrowCol[,'col']-1, '&y=', SrowCol[,'row']-1, '&z=',zoom,suffix, sep='')
        SrowCol$url = paste(path, SrowCol$tile, sep='')       
      } else if (tileNames == 'xyz') {
        SrowCol$cache = cacheDir2
        SrowCol$tile = paste(SrowCol[,'col']-1, '/', SrowCol[,'row']-1, '/',zoom,suffix, sep='')
        SrowCol$url = paste(path, SrowCol$tile, sep='')       
      } else if (tileNames == 'zyx') {
        SrowCol$cache = file.path(cacheDir2, zoom, SrowCol[,'row']-1)
        SrowCol$path = paste(path,zoom,'/',SrowCol[,'row']-1,'/',sep='')
        SrowCol$tile = paste(SrowCol[,'col']-1, suffix, sep='')
        SrowCol$url = paste(SrowCol$path, SrowCol$tile, sep='')
      } else {
        warning('tileNames must be zxy or zyx or xyz=')
      }
      
    SrowCol$file = file.path(SrowCol$cache , SrowCol$tile )
    SrowCol$index = 1:nrow(SrowCol)
  rasters = list()
#  colourtableList = list()
  
#  uniqueBlock = unique(SrowCol$block)
 # for(Dblock in uniqueBlock){
    
#    tilesHere = which(SrowCol$block == Dblock)
#    colourtableList[[Dblock]] = list()
#    rasters[[Dblock]] = list()

    for(Drow in SrowCol$index) {
      Dx = SrowCol[Drow, 'col']
      Dy = SrowCol[Drow, 'row']
      Dcell = SrowCol[Drow, 'cell']      
      Dcache = SrowCol[Drow, 'cache']
      Dfile = SrowCol[Drow, 'file']
      Durl = SrowCol[Drow, 'url']


      # extent of the cell
      Dextent = terra::ext(rep(terra::xyFromCell(rasterSphere, Dcell), each=2) + rep(terra::res(rasterSphere), each=2)*c(-1,1,-1,1)/2)
      
      

      dirCreateMapmisc(path=Dcache,recursive=TRUE,showWarnings=FALSE)
      
      
      Dsize = file.info(Dfile)['size']
      if(!any(Dsize > 0,na.rm=TRUE)) {
        if(verbose) cat("downloading ", Durl, "\n")
        try(downloadFileMapmisc(
            Durl, Dfile, quiet=!verbose, 
            method='auto', mode = 'wb'
          ), silent=TRUE)
      } else {
        if(verbose) cat("tile ", Dfile, " cached\n")
      }
      
      thisimage = try(suppressWarnings(terra::rast(Dfile)), silent=TRUE)
      
      if(any(class(thisimage)=='try-error')) {
        if(verbose) warning("tile ", Dfile, " cannot be loaded")
        
        thisimage = terra::rast(
          Dextent, nrows=256, ncols=256, crs=crsMerc
        )
        terra::values(thisimage) = NA
      } else {
        terra::crs(thisimage) = crsMerc
        terra::ext(thisimage) = Dextent
      }

      rasters[[Drow]] = thisimage
    } # end Drow

#    plot(worldMap);for(D in 1:length(rasters)) {plotRGB(rasters[[D]], add=TRUE, legend=FALSE)}      

    colourTableList = lapply(rasters, terra::coltab)
    for(D in 1:length(colourTableList)) {
      colourTableList[[D]] = as.data.frame(colourTableList[[D]][[1]])
      if(nrow(colourTableList[[D]])) colourTableList[[D]]$tile = D
    }
    colourtable = do.call(rbind, colourTableList)



    if(nrow(colourtable)) {
      theDup = duplicated(colourtable[,c('red','green','blue','alpha')])
      names(colourtable) = gsub("^value$", "value.old", names(colourtable))
      colourtable = cbind(colourtable, value = NA)
      colourtable[!theDup, 'value'] = seq(0, len=sum(!theDup), by=1)

      colourtableUnique = colourtable[!theDup, ]
      colourtableDup = colourtable[theDup, setdiff(names(colourtable), c('value'))]
      colourtableMerge = merge(colourtableDup, colourtableUnique[, setdiff(names(colourtableUnique), c('tile','value.old','row'))], 
        by = c('red','green','blue','alpha'), suffixes = c('.old','.new'))

      colourtableAll = rbind(colourtableUnique[names(colourtable)], colourtableMerge[,names(colourtable)])
      Stiles = sort(unique(colourtableAll$tile))

      colourtableFinal = colourtableUnique[order(colourtableUnique$value),c('value','red','green','blue','alpha')]

#      xx = rasters
      for(Dtile in Stiles) {
        toSub = colourtableAll[colourtableAll$tile == Dtile, ]
        rasters[[Dtile]] = terra::subst(rasters[[Dtile]], from=toSub[,'value.old'], to=toSub[,'value'])
        terra::coltab(rasters[[Dtile]]) = colourtableFinal
      }
    } # if have colourtable


    # merge into blocks

    SrowCol$newBlock = c(TRUE, !(diff(SrowCol$col)==0 & diff(SrowCol$row)==1))
    SrowCol$block = cumsum(SrowCol$newBlock)

    rastersMerged = list()
    for(Dblock in unique(SrowCol$block)) {
      blockHere = SrowCol[SrowCol$block == Dblock, ]
      toMerge = rasters[blockHere$index]
      if(length(toMerge) > 1) {
        names(toMerge)[1:2] = c('x', 'y')
        rastersMerged[[Dblock]] = do.call(terra::merge, toMerge)
      } else {
        rastersMerged[[Dblock]] = rasters[[blockHere$index]]
      }
    }
#    plot(worldMap);for(D in 1:length(rastersMerged)) {plotRGB(rastersMerged[[D]], add=TRUE)}      


    # fill in values of the cells for output raster
    if(verbose) cat('reprojecting\n')
    terra::nlyr(outraster) = terra::nlyr(rastersMerged[[1]])
    values(outraster) = NA
    xSeq = xFromCol(outraster)
    SrowColSub = SrowCol[,c('row','col','cell','index','block')]

    SoutRows = unique(c(seq(1, nrow(outraster), by=NrowsPerCycle), nrow(outraster)+1))
    SoutCells = cellFromRowCol(outraster, SoutRows, 1)
    SoutCells[length(SoutCells)] = ncell(outraster)+1

    if(verbose) cat('reprojecting:  ', length(SoutRows), ' cycles:')
    for(Dcycle in seq(1,length(SoutRows)-1) ) {
      if(verbose) cat(Dcycle, ' ')

      ScellOut = seq(SoutCells[Dcycle], SoutCells[Dcycle+1]-1)
      thisRow = terra::project(
        vect(xyFromCell(outraster, ScellOut), crs=crs(outraster)), 
        crsMerc)


      SrowColHere = cbind(
        indexOut = 1:length(ScellOut),
        ScellOut = ScellOut,
        cell = terra::extract(rasterSphere, thisRow, cells=TRUE)[,'cell'])

      SrowColHere = merge(SrowColHere, SrowColSub, all.x=TRUE, all.y=FALSE)
      SrowColHere = split(SrowColHere, SrowColHere$block)
      outValuesHere = lapply(SrowColHere, function(xx) {
          cbind(ScellOut = xx[,'ScellOut'], 
            terra::extract(
              rastersMerged[[xx[1,'block']]], 
              thisRow[xx[, 'indexOut']], cells=FALSE, xy=FALSE, ID=FALSE, raw=TRUE))
      })
      outValuesHere = do.call(rbind, outValuesHere)
      outValuesHere = outValuesHere[order(outValuesHere[,1]), ]
      if(inMemory(outraster)) {
        values(outraster)[outValuesHere[,1],] = outValuesHere[,-1]
      } else {
        writeValues(outraster, outValuesHere[,2], outValuesHere[1,1], 1)
      } 

    }
    if(verbose) cat(Dcycle, '\n')




  terra::coltab(outraster) = colourtableFinal


  if(terra::nlyr(outraster)== 3) {
    names(outraster) = c('red','green','blue')
    RGB(outraster) = 1:3
  }
  if(terra::nlyr(outraster)== 4) {
    names(outraster) = c('red','green','blue', 'alpha')
    RGB(outraster) = 1:4
  }

  attributes(outraster)$tiles = 
    list(tiles = nrow(SrowCol), 
      zoom=zoom,
      path=path)
  
  return(outraster)	
  
}


