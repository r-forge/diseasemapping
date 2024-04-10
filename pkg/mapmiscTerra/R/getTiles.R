
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
  
  worldNrcan = suppressWarnings(terra::project(worldLL, nrCrs, res=1000))
  
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
  
  terra::rast(ext(extentMerc), nrows = N, ncols=N, crs=crsMerc)
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
  
  SrowCol = getRowCol(extMerc, zoom=zoom)
  
  length(SrowCol[[1]])*length(SrowCol[[2]])
  
}


getTiles = function(
  outraster, 
  zoom=1, 
  path="http://tile.openstreetmap.org/",
  cachePath='.',
  cacheDir,
  verbose=FALSE, suffix = '.png',
  tileNames = 'zxy'){

  if(missing(cacheDir)) {
    cacheDir = make.names(gsub(
          "^http.*//([[:alpha:]][.])*((tile|basemap)s?[.][[:digit:]]?)?(openstreetmap[.])?|[[:punct:]]$", 
          "", path))
  }
  
  maxPixelsPerCycle = 1e5
  NtestCols = 100
  NrowsPerCycle = floor(maxPixelsPerCycle/terra::ncol(outraster)) # number of rows of out raster to process simultaneously

  cacheDir2 = file.path(cachePath, cacheDir)
  
  rasterSphere = .getRasterMerc(zoom)  
  
  samplePoints = rast(terra::ext(outraster), res= (terra::xmax(outraster)-terra::xmin(outraster))/NtestCols, crs=terra::crs(outraster))
  samplePoints = vect(terra::xyFromCell(samplePoints, 1:terra::ncell(samplePoints)), crs=terra::crs(outraster))
  xMerc = suppressWarnings(terra::project(samplePoints, crsMerc))


  SrowColFull = terra::cellFromXY(rasterSphere, terra::crds(xMerc))
#  SrowColFull = terra::extract(rasterSphere, xMerc, cells=TRUE)[,'cell']
  SrowColFull = cbind(cell = SrowColFull, 
    row = as.integer(terra::rowFromCell(rasterSphere, SrowColFull)), 
    col = as.integer(terra::colFromCell(rasterSphere, SrowColFull)))

  SrowColFull = SrowColFull[!is.na(SrowColFull[,1]), ,drop=FALSE]
  SrowCol = SrowColFull[!duplicated(SrowColFull[,c('row','col')]), ,drop=FALSE]#c('row','col')]
  SrowCol = SrowCol[order(SrowCol[,'col'], SrowCol[,'row']),,drop=FALSE]
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
    SrowCol$index = 1L:nrow(SrowCol)
    SrowCol$bad = NA
  rasters = list()
#  colourtableList = list()
  
#  uniqueBlock = unique(SrowCol$block)
 # for(Dblock in uniqueBlock){
    
#    tilesHere = which(SrowCol$block == Dblock)
#    colourtableList[[Dblock]] = list()
#    rasters[[Dblock]] = list()

  if(nrow(SrowCol) > 50) {
    warning("number of map tiles is large, zoom ", zoom, " and ", nrow(SrowCol), " tiles")
  }

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
        SrowCol[Drow, 'bad'] = TRUE
      } else {
        SrowCol[Drow, 'bad'] = FALSE
        terra::crs(thisimage) = crsMerc
        terra::ext(thisimage) = Dextent
      }

      rasters[[Drow]] = thisimage
    } # end Drow

    Snlyrs = unlist(lapply(rasters, terra::nlyr))
    if(any(Snlyrs==1) & any(Snlyrs > 1) ){
      # mixture of colour tables and rgb
      Srow = which(Snlyrs == 1)
      for(Drow in Srow) rasters[[Drow]] = terra::colorize(rasters[[Drow]], 'rgb')
    }

#    map.new(rasters[[1]], buffer=2*(xmax(rasters[[1]])-xmin(rasters[[1]])));plot(worldMap, add=TRUE);for(D in 1:length(rasters)) {plot(rasters[[D]], add=TRUE, legend=FALSE)}      

#    map.new(worldMap);plot(worldMap, add=TRUE);for(D in 1:length(rasters)) {plot(rasters[[D]], add=TRUE, legend=FALSE)}      

    colourTableList = lapply(rasters, terra::coltab)
    for(D in seq(1, len=length(colourTableList)) ) {
      colourTableList[[D]] = as.data.frame(colourTableList[[D]][[1]])
      if(nrow(colourTableList[[D]])) colourTableList[[D]]$tile = D
    }
    colourtable = do.call(rbind, colourTableList)



    if(nrow(colourtable)) {
      theDup = duplicated(colourtable[,c('red','green','blue','alpha')])
      names(colourtable) = gsub("^value$", "value.old", names(colourtable))
      colourtable = cbind(colourtable, value = as.integer(NA))
      colourtable[!theDup, 'value'] = seq(0L, len=as.integer(sum(!theDup)), by=1L)

      colourtableUnique = colourtable[!theDup, ,drop=FALSE]
      colourtableDup = colourtable[theDup, setdiff(names(colourtable), c('value')),drop=FALSE]
      colourtableMerge = merge(colourtableDup, colourtableUnique[, setdiff(names(colourtableUnique), c('tile','value.old','row')),drop=FALSE], 
        by = c('red','green','blue','alpha'), suffixes = c('.old','.new'))

      colourtableAll = rbind(colourtableUnique[names(colourtable)], colourtableMerge[,names(colourtable)])
      Stiles = sort(unique(colourtableAll$tile))

      colourtableFinal = colourtableUnique[order(colourtableUnique$value),c('value','red','green','blue','alpha')]

#      xx = rasters
      for(Dtile in Stiles) {
        toSub = colourtableAll[colourtableAll$tile == Dtile, ,drop=FALSE]

        terra::values(rasters[[Dtile]])[,1] = toSub[
          match(terra::values(rasters[[Dtile]])[,1], toSub[,'value.old']),
          'value']

# code below doesn't preserve integers
#        rasters[[Dtile]] = terra::subst(rasters[[Dtile]], 
#          from=toSub[,'value.old'], to=as.integer(toSub[,'value']), raw=TRUE)
        terra::coltab(rasters[[Dtile]]) = colourtableFinal
      }
    }  else {# if have colourtable
      colourtableFinal = NULL
    }
    # merge columns into blocks

    SrowCol$newBlock = c(TRUE, !(diff(SrowCol$col)==0 & diff(SrowCol$row)==1))
    SrowCol$block = as.integer(cumsum(SrowCol$newBlock))

    rastersMerged = list()
    for(Dblock in unique(SrowCol$block)) {
      blockHere = SrowCol[SrowCol$block == Dblock, ,drop=FALSE]
      toMerge = rasters[blockHere$index]
      if(length(toMerge) > 1) {
        names(toMerge)[1:2] = c('x', 'y')
        rastersMerged[[Dblock]] = do.call(terra::merge, toMerge)
      } else {
        rastersMerged[[Dblock]] = rasters[[blockHere$index]]
      }
#      plot(rastersMerged[[Dblock]], add=TRUE)
    }
#   map.new(worldMap);for(D in 1:length(rastersMerged)) {plot(rastersMerged[[D]], add=TRUE)};plot(worldMap,add=TRUE)      
    # merge rows if possible
    Nblocks = length(rastersMerged)
    if(Nblocks > 1) {
      toMerge = matrix(FALSE, Nblocks, Nblocks)
      for(DblockCol in seq(from=1,by=1, length=Nblocks-1)) { # won't do if Nblocks=1
        for(DblockRow in seq(DblockCol+1, Nblocks)) {
          bmat1 = SrowCol[SrowCol$block == DblockCol, ]
          bmat2 = SrowCol[SrowCol$block == DblockRow, ]
         if(nrow(bmat1) == nrow(bmat2)) {
            toMerge[DblockRow, DblockCol] = 
              all(bmat1$row == bmat2$row) & all(abs(bmat1$col - bmat2$col) ==1)
          }
        }
      }

   
    toMerge = toMerge + t(toMerge)
    toMerge = toMerge >= 1
    blockToMerge = blockToMergeOrig = which.max(apply(toMerge, 2, sum))
    for(D in seq(1, nrow(toMerge)^2)) {
      blockToMerge = unique(c(blockToMerge, unlist(apply(toMerge[, blockToMerge,drop=FALSE], 2, which))))
    }
    blockToMerge = sort(blockToMerge)
    blockToMergeOrig = min(blockToMerge)

    if(length(blockToMerge)>1) {
      rastersMerged[[blockToMergeOrig]] = do.call(terra::merge, rastersMerged[blockToMerge])
    }
#    rastersMerged = rastersMerged[sort(unique(c(blockToMergeOrig, setdiff(1:length(rastersMerged), blockToMerge))))]
    SrowCol[SrowCol$block %in% blockToMerge, 'block'] = blockToMergeOrig
    } # more than one block


    if(identical(crs(rastersMerged[[1]]), crs(outraster)) & (length(rastersMerged)==1) ) {
        outraster = rastersMerged[[1]]
    } else { # need to reproject


    # fill in values of the cells for output raster
    if(verbose) cat('reprojecting\n')
    terra::nlyr(outraster) = terra::nlyr(rastersMerged[[1]])
    terra::values(outraster) = as.integer(NA)
    xSeq = terra::xFromCol(outraster)
    SrowColSub = SrowCol[,c('row','col','cell','index','block'),drop=FALSE]

    SoutRows = unique(c(seq(1, terra::nrow(outraster), by=NrowsPerCycle), terra::nrow(outraster)+1))
    SoutCells = as.integer(terra::cellFromRowCol(outraster, SoutRows, 1))
    SoutCells[length(SoutCells)] = as.integer(terra::ncell(outraster))+1L

    if(verbose) cat('reprojecting:  ', length(SoutRows), ' cycles:')

#    terra::values(rasterSphere) = 1:terra::ncell(rasterSphere)
    for(Dcycle in seq(1,length(SoutRows)-1) ) {
      if(verbose) cat(Dcycle, ' ')

      ScellOut =  seq(SoutCells[Dcycle], SoutCells[Dcycle+1]-1) 
      thisRow = suppressWarnings(terra::project(
        vect(terra::xyFromCell(outraster, ScellOut), crs=crs(outraster)), 
        crsMerc, partial=TRUE))
      thisRowGeom = terra::geom(thisRow)[,c('x','y'), drop=FALSE]

      theCell = terra::cellFromXY(rasterSphere, thisRowGeom)

      SrowColHere = data.frame(
        indexOut = 1L:length(ScellOut),
        ScellOut = ScellOut,
        cell = theCell)

      SrowColHere = merge(SrowColHere, SrowColSub, all.x=TRUE, all.y=FALSE)
      SrowColHere = split(SrowColHere, SrowColHere$block)
      outValuesHere = lapply(SrowColHere, function(xx) {
            cellsHere = terra::cellFromXY(
              rastersMerged[[xx[1,'block']]], 
              crds(thisRow[xx[, 'indexOut',drop=FALSE]]))
            cbind(
              ScellOut = xx[,'ScellOut', drop=FALSE],
              terra::values(rastersMerged[[xx[1,'block']]])[cellsHere, , drop=FALSE])             
# doesn't preserve integers
#              values=terra::extract(
#                rastersMerged[[xx[1,'block']]], 
#                thisRow[xx[, 'indexOut',drop=FALSE]], cells=FALSE, xy=FALSE, ID=FALSE, raw=TRUE))
      })

      if(length(outValuesHere)) {
        for(D in names(outValuesHere)[-1]) {
          names(outValuesHere[[D]]) = names(outValuesHere[[1]])
        }
        outValuesHere = as.matrix(do.call(rbind, outValuesHere))
        outValuesHere = outValuesHere[order(outValuesHere[,1]), ]
        if(terra::inMemory(outraster)) {
          terra::values(outraster)[outValuesHere[,1],] = outValuesHere[,-1]
        } else {
          terra::writeValues(outraster, outValuesHere[,2], outValuesHere[1,1], 1)
        } 
      }
    }
    if(verbose) cat(Dcycle, '\n')

  } # end reproject

  
  terra::coltab(outraster) = colourtableFinal


  if(terra::nlyr(outraster)== 3) {
    names(outraster) = c('red','green','blue')
    terra::RGB(outraster) = 1:3
  }
  if(terra::nlyr(outraster)== 4) {
    names(outraster) = c('red','green','blue', 'alpha')
    terra::RGB(outraster) = 1:4
  }

  attributes(outraster)$tiles = 
    list(tiles = nrow(SrowCol), 
      zoom=zoom,
      path=path)
  
  return(outraster)	
  
}


