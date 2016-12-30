

crsModis <- CRS(
"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
	)

getModisRaster = function() {
modisRaster = raster(
		extent(-20015109.354,20015109.354,-10007554.677,10007554.677),
		nrow=18, ncol=36,
		crs=crsModis		
		)
		
values(modisRaster)= 1:ncell(modisRaster)
		modisDf = data.frame(
		ID = values(modisRaster),
		h=colFromCell(modisRaster, 1:ncell(modisRaster))-1,
		v=rowFromCell(modisRaster, 1:ncell(modisRaster))-1
		)
modisDf$hString = sprintf("%02d", modisDf$h)
modisDf$vString = sprintf("%02d", modisDf$v)

modisDf$tile = paste(
		'h', modisDf$hString, 'v', modisDf$vString, sep=''
		)

modisRaster = as.factor(modisRaster)
levels(modisRaster) = list(modisDf)
modisRaster
}

modisRaster = getModisRaster()

