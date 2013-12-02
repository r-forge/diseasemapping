library("geostatsp")
data("swissRain")
swissRain$elevation = extract(swissAltitude, swissRain)
swissRain$sqrtrain = sqrt(swissRain$rain)

# estimate parameters
swissFit =  likfitLgm(swissRain, trend=sqrtrain ~ elevation, 
		param=c(range=51700, nugget=0.11,shape=1,  
				anisoAngleDegrees=37, anisoRatio=7.6),
		paramToEstimate = c("range","nugget", 
				"anisoAngleDegrees", "anisoRatio"))

# simulate from the random effect conditional on
#   the observed data
swissSim = grfConditional(data=swissRain, y=swissFit$resid,
		param=swissFit$param, locations=30, 
		Nsim=4)

pdf("grfConditionalL.pdf",height=6,width=9)
par(mfrow=c(2,2),mar=c(0,0,0,0))
for(D in 1:4){
	plot(swissSim[[D]])
	plot(swissBorder,add=TRUE)
}
dev.off()



# create a small raster of elevation data
swissAltSmall = raster::aggregate(swissAltitude,fact=5)
# calculate the fixed effects portion of the rainfall process
rainMean = swissFit$param["(Intercept)"] +
			swissFit$param["elevation"] * swissAltSmall
rainMean = raster::mask(rainMean, swissBorder)
	
# define a function to identify the location of maximum rainfall	
maxRainLocation = function(x) {
	rain =  rainMean + x
	xyFromCell(rain, which.max(rain))
}



# get a conditional sample of three locations of maximum rainfall
swissRain$resid = swissFit$resid
swissLocation = grfConditional(data=swissRain, 
		y="resid",
		param=swissFit$param, locations=swissAltSmall, 
		Nsim=12, fun = maxRainLocation)


# plot the simulated random effect
pdf("grfConditionalL.pdf")
par(mar=c(0,0,0,0), mfrow=c(1,1))
plot(raster::overlay(swissSim, fun=mean)) 
plot(swissBorder, add=TRUE)
# add the locations to the map
points(t(swissLocation),col='red')
dev.off()


# more simulations
swissSim = grfConditional(data=swissRain, y=swissFit$resid,
		param=list(range = seq(20000, 50000, by=10000),
				variance=1, shape=1), locations=50, 
		Nsim=4)

pdf("grfConditionalParams.pdf",height=6,width=9)
par(mfrow=c(2,2),mar=c(0,0,0,0))
for(D in 1:4){
	plot(swissSim[[D]])
	plot(swissBorder,add=TRUE)
}
dev.off()

