# rongelap island

data(rongelapUTM)
rongelapLL = spTransform(rongelapUTM,  CRS("+init=epsg:4326"))

# Surface Reflectance 8-Day L3 Global 250m  http://e4ftl01.cr.usgs.gov/MOLA/
ndvi = raster('HDF4_EOS:EOS_GRID:\"/store/patrick/spatialData/rongelap/MYD09Q1.A2002185.h34v07.005.2007173043331.hdf\":MOD_Grid_250m_Surface_Reflectance:sur_refl_b02')
rongelapMODIS = spTransform(rongelapUTM,  ndvi@crs)

ndvi2 = crop(ndvi,extent(rongelapMODIS@bbox))

plot(ndvi2 )
plot(rongelapMODIS,add=T)

ndviUTM = projectRaster(ndvi2, crs= (rongelapUTM@proj4string))
names(ndviUTM) = "ndvi"
plot(ndviUTM)
points(rongelapUTM)
rongelapUTM$logtime = log(rongelapUTM$time)
res = glgm(data=rongelapUTM, formula = count ~ offset(logtime) + ndvi, covariates=ndviUTM, family="poisson",cells=48)
# data=rongelapUTM; formula = count ~ offset(logtime) + ndvi; covariates=ndviUTM; family="poisson";cells=30;priorCI=NULL;maternRoughness=1; buffer=0
 plot(res$raster[["predict.exp"]])