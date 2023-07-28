crsLL = terra::crs("+init=epsg:4326")
crsMerc = terra::crs("+init=epsg:3857")
crsModis <- terra::crs("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

crsCanada = terra::crs("+proj=omerc +lat_0=54.766 +lonc=-101.876 +alpha=-83.5 +k=0.998 +x_0=0 +y_0=0 +gamma=-83.498 +ellps=WGS84 +units=m")

bboxLL = terra::wrap(terra::vect(as.matrix(expand.grid(c(-180,180), c(-90,90))),
  type="polygons", crs=crsLL))

eps = 0.5

bboxLLsafe = terra::wrap(terra::vect(cbind(
    c(-180+eps,180-eps,180-eps, -180+eps), c(-90+eps,-90+eps,90-eps,90-eps)),
  type="polygons", crs=crsLL))
rm(eps)


extentMerc = c(xim = -20037508.3427892, xmax=20037508.3427892, ymin=-20037508.3427893, ymax=20037508.3427892)





