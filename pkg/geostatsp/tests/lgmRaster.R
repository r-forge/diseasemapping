library('geostatsp')
data("swissRainR")

anotherx = raster(swissRainR[['alt']])
values(anotherx) = seq(0,1,len=ncell(anotherx))
names(anotherx) = "myvar"

swissRainR2 = brick(swissRainR[['alt']], 
		sqrt(swissRainR[['prec1']]),
		anotherx)

swissResR =  lgm( formula=layer ~ alt+ myvar, 
		data=swissRainR2, shape=2,
		oneminusar=seq(0.05, 0.1, len=6),
		nugget =  seq(0,0.01,len=20),
		adjustEdges=FALSE,
		mc.cores=c(1,2)[1+(.Platform$OS.type=='unix')]
)

swissResR$summary[c('oneminusar','range','propNugget',
				grep("\\.betaHat$", rownames(swissResR$summary), value=TRUE)),]

# profile likelihood plot

# plot map of predicted values

swissResR =  lgm( formula=layer ~ alt+ myvar,  
		data=swissRainR2, shape=2,
		oneminusar=seq(0.05, 0.1, len=3),
		nugget =  seq(0,0.01,len=5),
		adjustEdges=TRUE,
		mc.cores=c(1,2)[1+(.Platform$OS.type=='unix')]
)

swissResR$summary[c('oneminusar','range','propNugget',
				grep("\\.betaHat$", rownames(swissResR$summary), value=TRUE)),]

# range in km
swissResR$summary[ 'range' ,] * sqrt(mean(values(area(swissRainR))))/mean(res(swissRainR))

# profile likelihood plot

# plot map of predicted values
