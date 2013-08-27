data=swissRain; trend=sqrtrain ~ elevation 
coordinates=data;upper=lower=parscale=NULL;
param=c(range=51000, nugget=0.1,rough=1,  
		aniso.angle.degrees=35, aniso.ratio=7)
paramToEstimate = c("range","nugget", 
		"aniso.angle.degrees", "aniso.ratio")
reml=T

trend = myTrend2
		param=myParams2
		locations = swissRaster
		covariates=NULL; nugget.in.prediction=FALSE

 param=myParams
trend=myTrend
locations=swissRaster; Nsim=2
covariates=swissAltitude; nugget.in.prediction=TRUE
fun=NULL

source("/home/patrick/workspace/diseasemapping/pkg/geostatsp/R/formulaLhs.R")