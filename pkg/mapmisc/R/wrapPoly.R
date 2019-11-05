wrapPoly = function(x, crs) {
  if (is.null(attributes(crs)$crop)) {
    attributes(crs)$crop = llCropBox(crs)$crop
  }
  
  if (requireNamespace('rgeos', quietly = TRUE) &
      requireNamespace('rgdal', quietly = TRUE)) {
    toCropX = spTransform(attributes(crs)$crop, crs(x))
    xCrop = rgeos::gDifference(x, toCropX, byid = TRUE)
    
    row.names(xCrop) = gsub(" (buffer|[[:digit:]]+)$", "", row.names(xCrop))
    
    
    # remove short line segments
    if (length(grep("^SpatialLines", class(xCrop)))) {
      NsegPerLine = unlist(lapply(xCrop@lines,
                                  function(qq)
                                    length(qq@Lines)))
      
      for (D in which(NsegPerLine > 3)) {
        # retain only three biggest segments
        lineD = xCrop@lines[[D]]@Lines
        NpointsPerSeg = unlist(lapply(lineD,
                                      function(qq)
                                        nrow(qq@coords)))
        lineD = lineD[which(NpointsPerSeg > 2)]
        xCrop@lines[[D]]@Lines = lineD
      }
    }
    
    if (any(slotNames(x) == 'data')) {
      xCropData = x@data[match(row.names(xCrop),
                               rownames(x@data)), ,
          drop = FALSE]
      
      rownames(xCropData) = names(xCrop)
      
      xCrop = SpatialPolygonsDataFrame(xCrop,
                                       data = xCropData)
      
    }
    
    xTcrop = spTransform(xCrop, crs)
    
  } else {
    xTcrop = NULL
  }
  
  xTcrop
}

llCropBox = function(crs, res = 1) {
  
   
  
  extraSeq = seq(-0.75,0.75, by=0.1)
  extraBox = expand.grid(apply(
    expand.grid(c(-179,179),extraSeq),
              1, sum),
    apply(
      expand.grid(c(-89,89),extraSeq),
            1, sum))
  
  

  llPoints = rbind(polyhedron@coords,
                   as.matrix(extraBox), llBorder@coords)
  
  if (!requireNamespace('rgdal', quietly = TRUE)) {
    warning("rgdal package is required for this operation")
    return(NULL)
  }
  
  pointsTransIn = suppressWarnings(rgdal::rawTransform(
    as.character(crsLL),
    as.character(crs),
    nrow(llPoints),
    llPoints[, 1],
    llPoints[, 2]
  ))
  
  pointsInRegion = is.finite(pointsTransIn[[1]]) &
    is.finite(pointsTransIn[[2]])
  
  transInRegion = cbind(pointsTransIn[[1]][pointsInRegion],
                        pointsTransIn[[2]][pointsInRegion])
  transInRegion = transInRegion[order(transInRegion[, 1], transInRegion[, 2]), ]
  
  # if omerc, truncate the x range
  if (length(grep("proj=omerc", as.character(crs)))) {
    #    crsUp = gsub("(alpha|gamma)=([[:digit:]]|[.])+", "", as.character(crs))
    crsUp = gsub("gamma=([[:digit:]]|[.])+", "gamma=0", 
                 as.character(crs))    
    
    transInRegionUp = suppressWarnings(
      rgdal::rawTransform(
        as.character(crs), 
        as.character(crsUp), nrow(transInRegion), 
        transInRegion[,1], transInRegion[,2]))
    toTrunc = transInRegionUp[[1]]
    toTrunc = toTrunc[is.finite(toTrunc)]
    toTrunc = quantile(toTrunc, prob = c(0.01, 0.99), na.rm=TRUE)
    getRid =  which(
      transInRegionUp[[1]] < toTrunc[1] |
        transInRegionUp[[1]] > toTrunc[2])
    if(length(getRid))
      transInRegion = transInRegion[-getRid, ]
  }
  transInRegion = SpatialPoints(transInRegion,
                                proj4string = crs)
  
  # region in crs
  regionTransOrig = 
    rgeos::gConvexHull(transInRegion, 
                       byid = FALSE)
  # border in crs
  border2 = SpatialLines(list(Lines(list(Line(
    regionTransOrig@polygons[[1]]@Polygons[[1]]@coords
    )), ID=1)), proj4string = crs)
  border3=rgeos::gInterpolate(border2,seq(0,1,len=4000),
                              normalized=TRUE)
  
  # border of crs transformed to LL
  borderLL = suppressWarnings(rgdal::rawTransform(
    as.character(crs),
    as.character(crsLL),
    nrow(border3@coords),
    border3@coords[, 1],
    border3@coords[, 2]
  ))
  
  borderLL = cbind(borderLL[[1]],borderLL[[2]])
  borderLL = borderLL[is.finite(borderLL[,1]), ]
  
  theBreaks = c(0,which(abs(diff(borderLL[,1]))>350), nrow(borderLL))
  theLines = list()
  for(D in 1:(length(theBreaks)-1)){
    theLines[[D]] = Line(borderLL[
      seq(theBreaks[D]+1, theBreaks[D+1]), ])
  }
  
  borderLL2 = SpatialLines(list(Lines(
    theLines, ID=1)), proj4string = crsLL)
  
  borderLL3=rgeos::gInterpolate(borderLL2,seq(0,1,len=4000),
                              normalized=TRUE)
  borderLL3@proj4string = CRS()
  holeLL = rgeos::gBuffer(borderLL3,
          width = res)
  holeLL@proj4string = crsLL
  
  # get rid of holes
  notHoles = which(!unlist(lapply(holeLL@polygons[[1]]@Polygons,
                           function(xx)
                             xx@hole)))
  edgeLL = SpatialPolygons(
    list(Polygons(holeLL@polygons[[1]]@Polygons[notHoles], ID=1))
  )
  

  llBorderT = suppressWarnings(rgdal::rawTransform(
    as.character(crsLL),
    as.character(crs),
    nrow(llBorder@coords),
    llBorder@coords[, 1],
    llBorder@coords[, 2]
  ))
  llBorderT = cbind(llBorderT[[1]], llBorderT[[2]])  
  llBorderT = llBorderT[is.finite(llBorderT[,1]), ]
  
  resTrans = res * mean(apply(bbox(regionTransOrig), 1, diff) * (0.25 /
                                                                   180))
  
  borderTrans = rgeos::gSimplify(rgeos::gBuffer(
    SpatialPoints(llBorderT), width =
                                                  4 * resTrans),
                                 tol = 4 * resTrans)
  projection(borderTrans) = crs
  
  regionTrans = rgeos::gSimplify(regionTransOrig,tol=resTrans)
  
    regionTransSmallInclude = rgeos::gDifference(
      regionTrans,
    borderTrans
  )
  
    # projectable region in LL
    transInRegion2 = rgeos::gIntersection(
      transInRegion, regionTransSmallInclude
    )
    pointsInLL = suppressWarnings(rgdal::rawTransform(
      as.character(crs),
      as.character(crsLL),
      nrow(transInRegion2@coords),
      transInRegion2@coords[, 1],
      transInRegion2@coords[, 2]
    ))
    pointsInLL2 = SpatialPoints(cbind(
      pointsInLL[[1]], pointsInLL[[2]]))
    
    # region in crs
    regionLLOrig = 
      rgeos::gConvexHull(pointsInLL2, 
                         byid = FALSE)
    regionLL = rgeos::gDifference(regionLLOrig, edgeLL)
    
    regionLL@proj4string = edgeLL@proj4string = crsLL
    
    
    
  result = list(
    crop = edgeLL,
    poly = regionLL,
    ellipse = regionTrans,
    polyTrans = regionTransSmallInclude
  )

}

