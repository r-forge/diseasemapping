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
  N = 20
  res = 1
  edge = c(0.05,1)
  oneSeq = seq(-1,1,len=N+1)[-1]-edge[1]/N
  coordsBox = expand.grid(
    x=oneSeq,y=oneSeq,z=oneSeq
    )
  coordsBox[,'r'] =  sqrt(apply(coordsBox^2,1,sum))
  coordsBox = coordsBox[coordsBox[,'r']<1,]
  coordsBox = coordsBox[coordsBox[,'r']>0,]
  coordsBox[,'lon'] = atan(coordsBox[,'y']/coordsBox[,'x'])
  coordsBox[,'lat'] = acos(coordsBox[,'z']/coordsBox[,'r'])/2-pi/4
  
  
  llPoints = SpatialPoints(
    as.matrix(coordsBox[,c('lon','lat')])*(360/(pi)),
    proj4string=crsLL
    )
  
  
    latSeq = sort(unique(llPoints@coords[,2]))
    lonSeq = sort(unique(llPoints@coords[,1]))
    
  llBorder = rbind(
    expand.grid(
      long=lonSeq, 
      lat=c(1,-1)*90-edge[2]
    ),
    expand.grid(
      long=c(1,-1)*180-edge[1], 
      lat=latSeq
    )
  )
  
  
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
  
  resTrans = mean(apply(bbox(regionTransOrig), 1, diff)*(res/180))
  borderTrans = rgeos::gBuffer(SpatialPoints(transOnBorder), width=1.2*resTrans)
  regionTransSmall = rgeos::gBuffer(regionTransOrig, width=-resTrans/2)
  
  # crop out areas which are close to edges in LL
  regionTransSmallInclude = rgeos::gSimplify(
    rgeos::gDifference(regionTransSmall, borderTrans),
    tol=resTrans/4, topologyPreserve=FALSE)
  # convert to separate polygons
  regionTransSmallInclude = regionTransSmallInclude@polygons[[1]]@Polygons
  regionTransSmallInclude = SpatialPolygons(
    mapply(function(srl, ID) Polygons(list(srl), ID),
    srl=regionTransSmallInclude, ID=1:length(regionTransSmallInclude)
  ), proj4string = crs)
  regionTransSmallInclude = regionTransSmallInclude[
    order(rgeos::gArea(regionTransSmallInclude, byid=TRUE),decreasing=TRUE)
    ,]
  
  edgeTrans = sp::spsample(as(regionTransOrig, 'SpatialLines'),
    n=2*pi*90/min(res), type='regular')
  
  
  regionTransSmallInclude@proj4string = borderTrans@proj4string = 
    regionTransSmall@proj4string = edgeTrans@proj4string = crs
  
  edgeLL = spTransform(edgeTrans, crsLL)
  edgeLL = rgeos::gBuffer(SpatialPoints(edgeLL@coords), width=max(res))
  edgeLL@proj4string = crsLL
  
  regionLL = spTransform(regionTransSmallInclude, crsLL)
  
  list(crop=edgeLL, poly=regionLL, ellipse = regionTransSmall, 
    polyTrans = regionTransSmallInclude)
  
}

