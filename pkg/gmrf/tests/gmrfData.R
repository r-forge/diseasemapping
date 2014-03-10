

dataDir = '/store/patrick/spatialData'
#' distance from coast
#' http://worldgrids.org/doku.php

download.file(
		'http://worldgrids.org/lib/exe/fetch.php?media=dicgsh0a.tif.gz',
		file.path(dataDir,'coast.tif.gz'))
library("R.utils")
gunzip(file.path(dataDir,'coast.tif.gz'))
coast = raster(file.path(dataDir,'coast.tif'))

# population

download.file(
		'http://worldgrids.org/lib/exe/fetch.php?media=pdmgpw1a.tif.gz',
		file.path(dataDir,'pop.tif.gz'))
gunzip(file.path(dataDir,'pop.tif.gz'))
pop = raster(file.path(dataDir,'pop.tif'))

# slope
download.file(
		'http://worldgrids.org/lib/exe/fetch.php?media=slpsrt1a.tif.gz',
		file.path(dataDir,'slope.tif.gz'))
gunzip(file.path(dataDir,'slope.tif.gz'))
slope = raster(file.path(dataDir,'slope.tif'))

# rain march
download.file(
		'http://worldgrids.org/lib/exe/fetch.php?media=p03dwd1a.tif.gz',
		file.path(dataDir,'rainmarch.tif.gz'))
gunzip(file.path(dataDir,'rainmarch.tif.gz'))
rainmarch = raster(file.path(dataDir,'rainmarch.tif'))

bcFile = file.path(dataDir, "blackCarbon.nc") 
#'

#+ blackCarbonDownload, eval=reDownload
download.file(
		'ftp://goldsmr1.sci.gsfc.nasa.gov/data/s4pa//GoCART_monthly/G4P0_1MS_2D_cc_aot.006/2004/g4p0e006tld_cc_aot2d_1MAVG_200401.nc',
		bcFile
)
#'

#' 'apb' is human-produced carbon

#+ blackCarbonLoad,cache=TRUE
bcLayers = MODIS::getSds(bcFile)
bcWhichLayers = grep("abp",bcLayers$SDSnames)
bcAll = raster::stack(bcLayers$SDS4gdal[bcWhichLayers])
extent(bcAll) = extent(-181.25,178.75,-91,91)
projection(bcAll) = CRS("+init=epsg:4326")


library('gmrf')
data('swissRainR')

swissRainDf = as.data.frame(swissRainR)
swissRainDf$alt = pmax(swissRainDf$alt,1)
for(D in 1:12) {
swissRainDf$rainTrans =log(swissRainDf[[paste("prec",D,sep="")]])
hist(swissRainDf$rainTrans)
locator(1)
#if(FALSE){
	library('mgcv')
	plot(gam(rainTrans~s(log(alt)),data=swissRainDf))
	abline(v=log(c(150,500,1500)))
	locator(1)
	#}



Xmat = cbind(intercept=1,
		altLow=pmin(log(swissRainDf$alt)-log(150),0), 
	altHigh=pmax(log(swissRainDf$alt)-log(150),0),
	altVeryHigh = pmax(log(swissRainDf$alt)-log(1500),0))

swissArgs = list(
		Yvec=swissRainDf[,'rainTrans'],
		Xmat=Xmat,
		NN=swissNN,
		maternShape=2
)


res1 = mapply(loglikGmrf,ar=seq(0.99, 0.999999,len=100), 
		MoreArgs=swissArgs)

plot(res1['ar',],res1['m2logL',],type='l')
locator(1)
}
Snugget= seq(0, 0.01, len=12)
propNugget = seq(0,0.01,len=12)


res2 = loglikGmrfNugget(ar=0.97, propNugget=Snugget,
		Yvec=swissRainDf$rainInv,
		Xmat=Xmat, NN=swissNN, maternShape=2)



res1 = loglikGmrfNugget(ar=0.97, propNugget=Snugget,
		Yvec=log(swissRainDf[,'prec5']),
		Xmat=Xmat, NN=swissNN, maternShape=2)

plot(t(res1[c('propNugget','logL'),]),type='l')

res2 = loglikGmrfNugget(ar=0.98, 
		propNugget=Snugget,
		Yvec=log(swissRainDf[,'prec5']),
		Xmat=Xmat, NN=swissNN, maternShape=2)

lines(t(res2[c('propNugget','logL'),]),col='blue')
legend('topright',lty=1,col=c('black','blue'),
		legend=c(0.97, 0.98),
		title='ar')


swissArgs = list(
		propNugget = seq(0,0.01,len=12),
		Yvec=log(swissRainDf[,'prec5']),
		Xmat=Xmat,
		NN=swissNN,
		maternShape=2
)

library('parallel')
swissResN = mcmapply(loglikGmrfNugget, 
		ar=seq(0.95, 0.999,len=8),
		MoreArgs=swissArgs,				
		mc.cores=4,SIMPLIFY='array')


Sprob = c(1, 0.9999, 0.999, 0.99, 0.95, 0.8, 0.5,   0) 
Schisq = qchisq(Sprob,df=2)
thebreaks = max(swissResN['logL',,]) - Schisq

thebreaks[thebreaks==-Inf] = min(swissResN['logL',,])-1

thecol = rev(RColorBrewer::brewer.pal( length(Sprob)-1,"Spectral")) 

image(
		swissResN['propNugget',,1],
		swissResN['ar',1,],
		swissResN['logL',,],
		breaks = thebreaks,
		col=thecol
)
library('mapmisc')
legendBreaks('right',breaks=Sprob,col=thecol)

swissSummary = summaryGmrfFitNugget(swissResN, swissArgs)
swissSummary$summary

