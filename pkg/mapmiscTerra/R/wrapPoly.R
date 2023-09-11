
#  library('Rpolyhedra')
#  mypoly = getPolyhedron(source='dmccooey', polyhedron.name = 'geodesic icosahedron pattern 223 [26,0]')
#  xx = mypoly$getState()$getVertices()   
#  xxSph = pracma::cart2sph(as.matrix(xx))
#  isohedron = (180/pi)*xxSph[,1:2]
#  isohedron[,1] = isohedron[,1]-mean(isohedron[,1])/2
#  save(isohedron, file='pkg/mapmiscTerra/data/isohedron.RData', compress='xz')


polesLLPolyFun = function(buffer.width, leftright=TRUE) {
 polesLLPoly = vect(
  c(terra::buffer(vect(cbind(0,90), crs=crsLL), width=buffer.width, quadsegs=100)  ,
    terra::buffer(vect(cbind(0,-90), crs=crsLL), width=buffer.width, quadsegs=100))
) 
   #  edges of south pole, because of antartica

poleSeq = cbind(180, seq(-90, -86.5, len=21))
poleSeq = rbind(poleSeq, cbind(179, rev(poleSeq[,2])))
polesLLPoly = terra::union(
  terra::union(polesLLPoly, vect(poleSeq, type='polygons', crs=crsLL)),  
  vect(poleSeq %*% diag(c(-1,1)), type='polygons', crs=crsLL)
)
poleSeq[,2] = abs(poleSeq[,2])
polesLLPoly = terra::union(
  terra::union(polesLLPoly, vect(poleSeq, type='polygons', crs=crsLL)),  
  vect(poleSeq %*% diag(c(-1,1)), type='polygons', crs=crsLL)
)



 if(leftright) {
 sides = terra::crds(terra::as.points(terra::buffer(
  vect(cbind(0, seq(-89.9, 89.9,len=300)), type='lines', crs=crsLL), 
  width=buffer.width)))
 sides = sides[sides[,1]>0,]
 sides = rbind(c(-0.1, sides[1,2]), sides, c(-0.1, sides[nrow(sides),2]))
 sidesRight = cbind(180-sides[,1], sides[,2])
 sidesLeft = cbind(-180+sides[,1], sides[,2])
 sidesLLPoly = vect(c(
  vect(sidesRight, type='polygons', crs=crsLL),
  vect(sidesLeft, type='polygons', crs=crsLL)
))
 result =  terra::aggregate(vect(c(polesLLPoly, sidesLLPoly)))
} else {
  result =  terra::aggregate(polesLLPoly)
}
result 

}

wrapPoly = function(x, crs, buffer.width = 100*1000) {

  if (is.null(attributes(crs)$crop) ) {
    attributes(crs)$crop = llCropBox(crs, buffer.width=buffer.width)$crop
  }
  toCrop = attributes(crs)$crop 
#  plot(toCrop)

  toCropX = project(x,crs(toCrop))

  xCrop = terra::erase(toCropX, toCrop)

#  xCrop = terra::crop(xCrop, terra::unwrap(bboxLLsafe))
  # plot(xCrop)


  project(xCrop, crs)
}

llCropBox = function(crs, 
  buffer.width=50*1000, densify.interval = 25*1000, 
  crop.distance = 2.1e7, crop.poles = FALSE, crop.leftright=FALSE,
  remove.holes=TRUE, cycles = 2, ellipse = NULL) {


 if(is.null(ellipse)) {
  utils::data('isohedron')
  isohedron[,2] = pmin(pmax(-89.99, isohedron[,2]), 89.99)

  bboxLLsafe = terra::unwrap(bboxLLsafe)
  LLborderInner = list()
   for(D in c(1,2,3)) {
     xx  = terra::buffer(bboxLLsafe, width=-D*buffer.width)
     LLborderInner[[as.character(D)]] = terra::crds(terra::densify(xx, densify.interval))
   }
   LLborderInner[[length(LLborderInner)+1]] = terra::crds(terra::densify(bboxLLsafe, densify.interval))

 LLpointsFull = vect(rbind(isohedron, do.call(rbind, LLborderInner)),crs=crsLL)
 LLpoints = terra::deepcopy(LLpointsFull)


for(DprojIter in 1:cycles) {
 suppressWarnings(pointsTransIn <- 
  terra::crds(project(LLpoints,  crs, partial=FALSE)))


 pointsInRegion = is.finite(pointsTransIn[[1]]) &
  is.finite(pointsTransIn[[2]]) & abs(pointsTransIn[,1]) < crop.distance  &
    abs(pointsTransIn[,2]) < crop.distance 


 transInRegion =  pointsTransIn[pointsInRegion,]
 transInRegion = transInRegion[order(transInRegion[, 1], transInRegion[, 2]), ]
 transInRegion = vect(transInRegion, crs=crs)


  # region in crs
  regionTransOrig = terra::convHull(transInRegion)
  regionTransOrigCoords = terra::crds(regionTransOrig)
  if(nrow(regionTransOrigCoords) < 100){
    regionTransOrig = terra::densify(regionTransOrig, densify.interval)
    regionTransOrigCoords = terra::crds(regionTransOrig)
  }

 theCentroid = apply(apply(regionTransOrigCoords, 2, range),2,mean)
 regionTransCentred = regionTransOrigCoords - matrix(theCentroid, nrow(regionTransOrigCoords), 2, byrow=TRUE)

 regionTransPolar = regionTransCentred[,1] + 1i*regionTransCentred[,2]

 # smooth a second time with angles 0 to 2pi, to smooth out end points
 regionTransPolarAngle2 = Arg(regionTransPolar)
  regionTransPolarAngle2neg =  regionTransPolarAngle2 < 0
  regionTransPolarAngle2[regionTransPolarAngle2neg] =   regionTransPolarAngle2[regionTransPolarAngle2neg]+2*pi

# if(length(regionTransPolar)>20) {
  Nout = 2001
  newxInner = seq(-pi/2, pi/2, len=Nout)
  newxEnds = seq(pi/2, 3*pi/2, len=Nout)

  smoothPolar = stats::smooth.spline(
    Arg(regionTransPolar), 
    Mod(regionTransPolar), all.knots=TRUE, 
    df=ceiling(length(regionTransPolar)*0.5))

  smoothPolar2 = stats::smooth.spline(
    regionTransPolarAngle2, 
    Mod(regionTransPolar), all.knots=TRUE, 
    df=ceiling(length(regionTransPolar)*0.5))



  polarDense = 0.999*c(stats::predict(smoothPolar, newxInner)$y * exp(1i*newxInner), 
    stats::predict(smoothPolar2, newxEnds)$y * exp(1i*newxEnds))
#  }

 polarDense = polarDense[order(Arg(polarDense))] + theCentroid[1] + 1i*theCentroid[2]

 regionTransOrigPoints = terra::vect(cbind(Re(polarDense), Im(polarDense)), crs=crs, type='points')

  LLpoints1 = suppressWarnings(project(regionTransOrigPoints, crsLL))
  LLpoints = terra::as.points(terra::buffer(LLpoints1, buffer.width, quadsegs = 5))
}

regionTransSmooth = terra::vect(cbind(Re(polarDense), Im(polarDense)), crs=crs, type='polygons')
regionTransSmooth =terra::densify( regionTransSmooth, interval = densify.interval)


} else { # have ellipse
  regionTransSmooth = ellipse # terra::densify(ellipse, densify.interval)
}

regionTransPoly1 = terra::densify(
  terra::buffer(regionTransSmooth, - 0.2*buffer.width, quadsegs = 10), 
  interval=densify.interval)
regionTransPoly2 = terra::densify(
  terra::buffer(regionTransSmooth, - 0.5* buffer.width, quadsegs = 10), 
  interval=densify.interval)


  # border of crs transformed to LL
borderLL1 = terra::geom(suppressWarnings(project(terra::as.points(regionTransPoly1), crsLL)))
borderLL1 = terra::vect(borderLL1[!is.nan(borderLL1[,'x']), ], crs=crsLL)

borderLL2 = terra::geom(suppressWarnings(project(terra::as.points(regionTransPoly2), crsLL)))
borderLL2 = terra::vect(borderLL2[!is.nan(borderLL2[,'x']), ], crs=crsLL)



# data('worldMap');worldMap = unwrap(worldMap);plot(project(worldMap, crsLL), ylim = c(-95,95));points(borderLL1, cex=0.3,col='blue');points(borderLL2, cex=0.1, col='red')     

whereIsJump1a = which(abs(diff(terra::crds(borderLL1)[,1])) > 180)
whereIsJump2a = which(abs(diff(terra::crds(borderLL2)[,1])) > 180)
whereIsJump1b = which(abs(diff(terra::crds(borderLL1)[,2])) > 100)
whereIsJump2b = which(abs(diff(terra::crds(borderLL2)[,2])) > 100)

theBreaks1 = diff(sort(unique(c(0,whereIsJump1a,whereIsJump1b, length(borderLL1)))))
theBreaks2 = diff(sort(unique(c(0,whereIsJump2a,whereIsJump2b, length(borderLL2)))))

borderLLsplit1 = split(borderLL1, rep(1:length(theBreaks1), theBreaks1))
borderLLsplit2 = split(borderLL2, rep(1:length(theBreaks2), theBreaks2))

borderLLsplit1 = borderLLsplit1[unlist(lapply(borderLLsplit1, function(xx) nrow(terra::crds(xx))))>0]
borderLLsplit2 = borderLLsplit2[unlist(lapply(borderLLsplit2, function(xx) nrow(terra::crds(xx))))>0]

# data(worldMap);worldMap = unwrap(worldMap);plot(project(worldMap, crsLL), ylim = c(-92, 92));for(D in 1:length(borderLLsplit2)) {plot(as.lines(borderLLsplit2[[D]]), add=TRUE, col=1+D, lwd=3)}   

borderLLsplit1l = lapply(borderLLsplit1, terra::as.lines)
borderLLsplit2l = lapply(borderLLsplit2, terra::as.lines)

borderLLsplitL = terra::vect(terra::svc(c(borderLLsplit1l, borderLLsplit2l)))
borderLLsplitL = terra::simplifyGeom(borderLLsplitL, tolerance = 0.05, preserveTopology=FALSE)
borderLLlinesListDens = terra::densify(borderLLsplitL, densify.interval)

allPoints =terra::as.points(borderLLlinesListDens)

if(crop.poles) {
  polesAndSidesLLpoly=polesLLPolyFun(buffer.width, crop.leftright)
  safePoints = terra::erase(allPoints, polesAndSidesLLpoly)
  edgeLL = terra::aggregate(terra::buffer(safePoints, buffer.width))
  edgeLL = terra::aggregate(terra::union(edgeLL, terra::buffer(polesAndSidesLLpoly, 2*buffer.width) )) 

} else {  
  edgeLL = terra::aggregate(terra::buffer(allPoints, buffer.width))
}


  if(remove.holes) {
    eps = 0.1
      leftSide = terra::ext(terra::crop(allPoints, terra::ext(-180,-180+eps,-90,90)))
      rightSide = terra::ext(terra::crop(allPoints, terra::ext(180-eps,180,-90,90)))
      if(length(leftSide)) {
        leftSide = terra::densify(vect(cbind(-180, c(terra::ymin(leftSide), terra::ymax(leftSide))), crs=crsLL, type='lines'), densify.interval)
        edgeLL = terra::aggregate(terra::union(edgeLL, terra::buffer(leftSide, buffer.width)))
      }
      if(length(rightSide)) {
        rightSide = terra::densify(vect(cbind(180, c(terra::ymin(rightSide), terra::ymax(rightSide))), crs=crsLL, type='lines'), densify.interval)
        edgeLL = terra::aggregate(terra::union(edgeLL, terra::buffer(rightSide, buffer.width)))
      }


    edgeLL = terra::fillHoles(edgeLL)

  }

# plot(project(worldMap, crsLL), ylim = c(-92, 92));plot(edgeLL, add=TRUE, col='blue')


regionCrs = rast(regionTransSmooth, nrow=100, ncol=100)
regionCrs = terra::rasterize(regionTransSmooth, regionCrs)
regionCrs = vect(c(terra::as.points(regionTransSmooth), terra::as.points(regionCrs)))
regionLL = project(regionCrs, crsLL)


return(list(
  crop = edgeLL,
  ellipse = regionTransSmooth,
  regionLL = regionLL
))

}

