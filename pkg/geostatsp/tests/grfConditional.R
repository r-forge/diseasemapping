library("geostatsp")
data("swissRain")
swissRain$elevation = extract(swissAltitude, swissRain)
swissRain$sqrtrain = sqrt(swissRain$rain)

swissFit =  likfitLgm(swissRain, trend=sqrtrain ~ elevation, 
		param=c(range=51700, nugget=0.11,rough=1,  
				aniso.angle.degrees=37, aniso.ratio=7.6),
		paramToEstimate = c("range","nugget", 
				"aniso.angle.degrees", "aniso.ratio"))
myTrend = swissFit$trend
myParams = swissFit$param

swissRaster = raster(extent(swissBorder), ncols=20, nrows=20, 
		crs=swissRain@proj4string)	

swissSim = grfConditional(data=swissRain, param=myParams,
		trend=myTrend, 
		locations=swissRaster, Nsim=2, 
		covariates=swissAltitude, nugget.in.prediction=TRUE) 

plot(mask(swissSim[[2]], swissBorder)^2) 
plot(swissBorder, add=TRUE)

