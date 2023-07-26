crsLL = terra::crs("+init=epsg:4326")
crsMerc = terra::crs("+init=epsg:3857")

crsModis <- terra::crs("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

bboxLL = terra::wrap(vect(as.matrix(expand.grid(c(-180,180), c(-90,90))),
  type="polygons", crs=crsLL))

eps = 0.5

bboxLLsafe = terra::wrap(vect(cbind(
    c(-180+eps,180-eps,180-eps, -180+eps), c(-90+eps,-90+eps,90-eps,90-eps)),
  type="polygons", crs=crsLL))



eps = 0.5

llBorder = rbind(
    cbind(seq(-180+eps, 180-eps, len=100), -90+eps), 
    cbind(180-eps, seq(-90+2*eps, 90-2*eps, len=100) ),
    cbind(seq(180-eps, -180+eps, len=100), 90-eps),
    cbind(-180, seq(90-2*eps, -90+2*eps, len=100))
  )
rm(eps)
gc()




