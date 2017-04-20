wrapPoly = function(x, crs){
  
  if(is.null(attributes(crs)$crop)) {
    attributes(crs)$crop = llCropBox(crs)$crop
  }
  
  if(requireNamespace('rgeos', quietly=TRUE) & 
    requireNamespace('rgdal', quietly=TRUE)) {	
    toCropX = spTransform(attributes(crs)$crop, crs(x))
    xCrop = rgeos::gDifference(x, toCropX, byid=TRUE)
    
    row.names(xCrop) = gsub(" (buffer|[[:digit:]]+)$","", row.names(xCrop))
    
    
    # remove short line segments
    if(length(grep("^SpatialLines", class(xCrop)))){
      
      NsegPerLine = unlist(lapply(xCrop@lines, 
          function(qq) length(qq@Lines)))
      
      for(D in which(NsegPerLine > 3)) {
        # retain only three biggest segments
        lineD = xCrop@lines[[D]]@Lines
        NpointsPerSeg = unlist(lapply(lineD, 
            function(qq) nrow(qq@coords)))
        lineD = lineD[which(NpointsPerSeg > 2)]
        xCrop@lines[[D]]@Lines = lineD
      }
    }
    
    if(any(slotNames(x)=='data')) {
      
      xCropData = x@data[match(
          row.names(xCrop),
          rownames(x@data)
        ),]
      
      rownames(xCropData) = names(xCrop)
      
      xCrop = SpatialPolygonsDataFrame(
        xCrop,
        data=xCropData
      )
      
    }
    
    xTcrop = spTransform(xCrop, crs)
    
  } else {
    xTcrop = NULL
  }
  
  xTcrop
}

llCropBox = function(crs, res=1) {
  
  # long-lat grid covering the globe
  N = 51
  res = 1
  edge = c(0.05,1)
  
  
  eps = 1/N
  oneSeq = seq(-1+2*eps,1-eps,len=N)
  
  coordsBox = expand.grid(
    x=oneSeq,y=oneSeq,z=oneSeq
  )
  coordsBox = coordsBox + sin(seq(0-eps,2*pi, len=nrow(coordsBox)))/N  
  coordsBox = coordsBox + 
    stats::runif(length(coordsBox), -0.6/N, 0.6/N)
  coordsBox[,'r'] =  sqrt(apply(coordsBox^2,1,sum))
  coordsBox = coordsBox[coordsBox[,'r'] < 1, ]
  coordsBox = coordsBox[coordsBox[,'r']> 0,]
  coordsBox[,'lon'] = atan(coordsBox[,'y']/coordsBox[,'x'])
  coordsBox[,'lat'] = acos(coordsBox[,'z']/coordsBox[,'r'])/2-pi/4
  
  
  llPoints =
    as.matrix(na.omit(coordsBox[,c('lon','lat')]))*(360/(pi))
  
  Sprob = seq(0,1,len=251)
  latSeq = sort(unique(c(
        seq(-90,90,len=201),
        quantile(llPoints[,2], prob = Sprob))))
  lonSeq = seq(-180,180,len=51)
  
  llBorder = cbind(
    lon=c(
      lonSeq, rep(-180, length(latSeq)),
      rep(180, length(latSeq))
    ),
    lat=c(
      rep(-90,length(lonSeq)), latSeq, latSeq
    ))
  
  
  if(!requireNamespace('rgdal', quietly=TRUE)) {
    warning("rgdal package is required for this operation")
    return(NULL)
  }
  pointsTransIn = suppressWarnings(
    rgdal::rawTransform(
      as.character(crsLL),
      as.character(crs),
      nrow(llPoints),
      llPoints[,1], 
      llPoints[,2]))
  
  pointsTransBorder = suppressWarnings(
    rgdal::rawTransform(
      as.character(crsLL),
      as.character(crs),
      nrow(llBorder),
      llBorder[,1], 
      llBorder[,2]))
  
  pointsInRegion = is.finite(
      pointsTransIn[[1]]) &
    is.finite(pointsTransIn[[2]])
  
  transInRegion = cbind(
    pointsTransIn[[1]][pointsInRegion],
    pointsTransIn[[2]][pointsInRegion])
  transInRegion = transInRegion[order(transInRegion[,1], transInRegion[,2]),]
  
  # if omerc, truncate the x range
  if(length(grep("proj=omerc", as.character(crs)))) {
    
    absX = abs(transInRegion[,1])
    toTrunc = quantile(absX, prob=c(0.95))
    transInRegion = transInRegion[absX < toTrunc,]
  }
  
  transOnBorder = cbind(
    pointsTransBorder[[1]],
    pointsTransBorder[[2]]
  )
  transOnBorder = transOnBorder[
    apply(transOnBorder, 1, function(xx) all(is.finite(xx))),
  ]    
  
  regionTransOrig = rgeos::gConvexHull(
    rgeos::gConvexHull(SpatialPoints(transInRegion), byid=FALSE)
  )
  resTrans = mean(apply(bbox(regionTransOrig), 1, diff)*(0.25/180))
  regionTransSmall = rgeos::gBuffer(regionTransOrig, width=-2*resTrans)
  
  
  if(nrow(transOnBorder)) {
    borderTrans = rgeos::gSimplify(
      rgeos::gBuffer(SpatialPoints(transOnBorder), width=2*resTrans),
      tol = 2*resTrans)
    # crop out areas which are close to edges in LL
    regionTransSmallInclude = #rgeos::gSimplify(
      rgeos::gDifference(regionTransSmall, 
        borderTrans)#,
#      tol=resTrans/4, topologyPreserve=FALSE)
  } else {
    borderTrans = SpatialPolygons(list())
    regionTransSmallInclude = regionTransSmall
  }
  
  # convert to separate polygons
  anyHoles = unlist(lapply(regionTransSmallInclude@polygons[[1]]@Polygons,
      function(xx) xx@hole))
  if(!any(anyHoles)) {
    regionTransSmallInclude = regionTransSmallInclude@polygons[[1]]@Polygons
    regionTransSmallInclude = SpatialPolygons(
      mapply(function(srl, ID) Polygons(list(srl), ID),
        srl=regionTransSmallInclude, ID=1:length(regionTransSmallInclude)
      ), proj4string = crs)
    regionTransSmallInclude = regionTransSmallInclude[
      order(rgeos::gArea(regionTransSmallInclude, byid=TRUE),decreasing=TRUE)
      ,]
  }
  edgeTrans = sp::spsample(as(regionTransOrig, 'SpatialLines'),
    n=3000, type='regular')
  
  regionTransSmallInclude@proj4string = borderTrans@proj4string = 
    regionTransOrig@proj4string = regionTransSmall@proj4string = 
    edgeTrans@proj4string = crs
  
  edgeLLP = spTransform(edgeTrans, crsLL)
  
  edgeLL = rgeos::gBuffer(SpatialPoints(edgeLLP@coords), width=res)
#  edgeLL = rgeos::gSimplify(edgeLL, tol=0.5)
  edgeLL@proj4string = crsLL
  
  regionLL = spTransform(regionTransSmallInclude, crsLL)
#  regionLL = raster::crop(regionLL, extent(-179.9, 179.9, -89.9, 89.9))
  
  list(crop=edgeLL, poly=regionLL, ellipse = regionTransSmall, 
    polyTrans = regionTransSmallInclude)
  
}

