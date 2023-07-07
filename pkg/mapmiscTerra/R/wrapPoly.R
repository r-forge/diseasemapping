
if(FALSE) {
library(Rpolyhedra)
mypoly = getPolyhedron(source='dmccooey', polyhedron.name = 'geodesic icosahedron pattern 223 [26,0]')
xx = mypoly$getState()$getVertices()   
xxSph = pracma::cart2sph(as.matrix(xx))
isohedron = (180/pi)*xxSph[,1:2]
isohedron[,1] = isohedron[,1]-mean(isohedron[,1])/2
save(isohedron, file='pkg/mapmiscTerra/data/isohedron.RData', compress='xz')
}


polesLLPolyFun = function(buffer.width) {
   polesLLPoly = vect(
    c(buffer(vect(cbind(0,90), crs=crsLL), buffer.width, quadsegs=100)  ,
      buffer(vect(cbind(0,-90), crs=crsLL), buffer.width, quadsegs=100)))
  sides = terra::crds(terra::as.points(terra::buffer(vect(  cbind(0, seq(-89.9, 89.9,len=300)), type='lines', crs=crsLL), buffer.width)))
  sides = sides[sides[,1]>0,]
  sides = rbind(c(-0.1, sides[1,2]), sides, c(-0.1, sides[nrow(sides),2]))
  sidesRight = cbind(180-sides[,1], sides[,2])
  sidesLeft = cbind(-180+sides[,1], sides[,2])
  sidesLLPoly = vect(c(
    vect(sidesRight, type='polygons', crs=crsLL),
    vect(sidesLeft, type='polygons', crs=crsLL)
  ))

   terra::aggregate(vect(c(polesLLPoly, sidesLLPoly)))
}

wrapPoly = function(x, crs, buffer.width) {

  if (is.null(attributes(crs)$crop) | !missing(buffer.width)) {
    attributes(crs)$crop = llCropBox(crs, buffer.width=buffer.width)$crop
  }
  toCrop = attributes(crs)$crop 
#  plot(toCrop)

  toCropX = project(x,crs(toCrop))

  xCrop = terra::erase(toCropX, toCrop)

  xCrop = terra::crop(xCrop, bboxLLsafe)
  # plot(xCrop)


 project(xCrop, crs)
}

llCropBox = function(crs, buffer.width=100*1000, densify.interval = 25*1000) {
  
   
  polesAndSidesLLpoly=polesLLPolyFun(buffer.width)

  data('isohedron')
  isohedron[,2] = pmin(pmax(-89.9, isohedron[,2]), 89.9)

  LLborderInner = list()
  for(D in c(2,5,10, 15)) {
   xx  = terra::buffer(vect(llBorder, type='polygon', crs=crsLL), width=-D*buffer.width)
   LLborderInner[[as.character(D)]] = terra::crds(terra::densify(xx, densify.interval))
  }
  LLpoints = vect(rbind(isohedron, llBorder, do.call(rbind, LLborderInner)),crs=crsLL)
  
  suppressWarnings(pointsTransIn <- 
      terra::crds(project(LLpoints,  crs, partial=FALSE)))

  
  pointsInRegion = is.finite(pointsTransIn[[1]]) &
    is.finite(pointsTransIn[[2]]) & abs(pointsTransIn[,1]) < 2e7 &
       abs(pointsTransIn[,2]) < 2e7
  

  transInRegion =  pointsTransIn[pointsInRegion,]
  transInRegion = transInRegion[order(transInRegion[, 1], transInRegion[, 2]), ]
  transInRegion = vect(transInRegion, crs=crs)


  # region in crs
  regionTransOrig = terra::convHull(transInRegion)
  regionTransOrigPoints = terra::as.points(regionTransOrig)

  if(length(regionTransOrigPoints)>100) {
  Nout = 5000
  smoothX = stats::smooth.spline(terra::crds(regionTransOrigPoints)[,1], df=length(regionTransOrigPoints)*0.5)
  smoothY = stats::smooth.spline(terra::crds(regionTransOrigPoints)[,2], df=length(regionTransOrigPoints)*0.5)

  newx = seq(1, length(regionTransOrigPoints), len=Nout)
  regionTransSmooth = terra::densify(vect(
    cbind(stats::predict(smoothX, newx)$y, stats::predict(smoothY, newx)$y), 
    type='polygons', crs=crs), interval = densify.interval)
} else {
  regionTransSmooth = terra::densify(regionTransOrig, densify.interval)
}

  regionTransPoly1 = terra::densify(
      terra::buffer(regionTransSmooth, -buffer.width/8), 
    interval=densify.interval)
  regionTransPoly2 = terra::densify(
      terra::buffer(regionTransSmooth, -buffer.width/5), 
    interval=densify.interval)


  # border of crs transformed to LL
  borderLL1 = suppressWarnings(project(terra::as.points(regionTransPoly1), crsLL))
  borderLL2 = suppressWarnings(project(terra::as.points(regionTransPoly2), crsLL))
  
  


# plot(project(worldMap, crsLL), ylim = c(-95,95));points(borderLL1, cex=0.1,col='blue');points(borderLL2, cex=0.1, col='red')     

  whereIsJump1a = which(abs(diff(terra::crds(borderLL1)[,1])) > 180)
  whereIsJump2a = which(abs(diff(terra::crds(borderLL2)[,1])) > 180)
  whereIsJump1b = which(abs(diff(terra::crds(borderLL1)[,2])) > 100)
  whereIsJump2b = which(abs(diff(terra::crds(borderLL2)[,2])) > 100)

  theBreaks1 = diff(sort(unique(c(0,whereIsJump1a,whereIsJump1b, length(borderLL1)))))
  theBreaks2 = diff(sort(unique(c(0,whereIsJump2a,whereIsJump2b, length(borderLL2)))))

  borderLLsplit1 = split(borderLL1, rep(1:length(theBreaks1), theBreaks1))
  borderLLsplit2 = split(borderLL2, rep(1:length(theBreaks2), theBreaks2))

# plot(project(worldMap, crsLL), ylim = c(-92, 92));for(D in 1:length(borderLLsplit2)) {plot(as.lines(borderLLsplit2[[D]]), add=TRUE, col=1+D, lwd=3)}   


  borderLLlinesList = c(lapply(borderLLsplit1, terra::as.lines), lapply(borderLLsplit2, terra::as.lines))
  borderLLlinesListDens = lapply(borderLLlinesList, function(xx) {
       terra::densify(xx, densify.interval)
  })
  allPoints = vect(as.matrix(do.call(rbind, lapply(borderLLlinesListDens, crds))), crs=crsLL)

  pointsOnBorders = terra::crop(allPoints, polesAndSidesLLpoly)
  safePoints = terra::erase(allPoints, polesAndSidesLLpoly)
  edgeLL = terra::aggregate(terra::buffer(safePoints, buffer.width))
  

  if(length(pointsOnBorders)>10) {
    polesAndSidesLLpolySplit = terra::split(polesAndSidesLLpoly, 
      vect(c(vect(cbind(c(-180,180),c(-90,90)), type='lines',crs=crsLL),vect(cbind(c(-180,180),c(90,-90)), type='lines',crs=crsLL))))
    NperQuadrant = as.data.frame(table(terra::extract(polesAndSidesLLpolySplit,pointsOnBorders)[,2]))
    toKeep = as.numeric(NperQuadrant[NperQuadrant[,2] > 10 , 1])
    if(length(toKeep))
      edgeLL = terra::aggregate(vect(c(edgeLL, polesAndSidesLLpolySplit[toKeep])))
  }

# plot(project(worldMap, crsLL), ylim = c(-92, 92));plot(edgeLL, add=TRUE, col='blue')

  if(!all(pointsInRegion)) {
    excludedFromLL =  LLpoints[!pointsInRegion]
    exclBuf = terra::buffer(terra::erase(excludedFromLL, polesAndSidesLLpoly), buffer.width)
    theHullLL = terra::convHull(vect(c(excludedFromLL, terra::as.points(exclBuf))))
    edgeLL =  terra::aggregate(terra::vect(c(edgeLL, theHullLL)))
  } 



  return(list(
    crop = edgeLL,
    ellipse = terra::as.polygons(terra::as.lines(regionTransPoly2))
  ))

}

if(FALSE) {
  holeLL = rgeos::gBuffer(borderLL3,
          width = res)
  holeLL@proj4string = bboxLL@proj4string
  holeLL = rgeos::gIntersection(
    holeLL, bboxLL)
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
  crs(borderTrans) = crs
  
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
    pointsInLL2 = vect(cbind(
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

