#' # LGCP

#+ toCompile, eval=FALSE, include=FALSE
knitr::spin("lgcpTest.R", format='Rmd', knit=TRUE)
#'

#+ packages
library('geostatsp')
library('diseasemapping')
library("mapmisc")
library("knitr")
if(file.exists("../docs/knitrCode.R")) {
	source("../docs/knitrCode.R")
	knit_hooks$set(plot=hook_plot_p) 
	opts_chunk$set(fig.ncol=2, fig.height=4, fig.width=8)
}
#'


#+ kentucky, cache=TRUE
data('kentucky')
data('kentuckyTract')
kentucky = spTransform(kentucky, mapmisc::omerc(kentucky))
kentuckyT = spTransform(kentuckyTract, projection(kentucky))
kMap = openmap(kentuckyT, path='stamen-toner')
kMap = tonerToTrans(kMap, power=1)
#'

#+ params
output.ras = squareRaster(kentuckyT, 100)
cov.ras = list()
cov.ras$w1 = cov.ras$w2 = output.ras
values(cov.ras$w1) = yFromCell(output.ras, 1:ncell(output.ras))/500000
values(cov.ras$w2) = xFromCell(output.ras, 1:ncell(output.ras))/500000
#'


#+ offsets, cache=TRUE
larynxRates= cancerRates("USA", year=1998:2002,site="Larynx")
larynxRates = larynxRates * 5
kentuckyT = getSMR(kentuckyT, larynxRates, regionCode="County")
kentuckyOffset = rasterize(
		kentuckyT,
		output.ras,
		field='logExpected_surfaceArea'
)
cov.ras$offset = kentuckyOffset 
#'

#+ simB, cache=TRUE, stuff=2
U = c(mean = 0, 
		variance=0.25^2, 
		range=120000,
		shape=2)

set.seed(0)
lgcp.sim = simLgcp(U,cov.ras,betas=c(-0.1,0.025),
		rasterTemplate=output.ras, n=4,
		offset='offset')
#'

#+ forSimPlot
Ssim = grep('^linear|^intensity|^offset|^w[[:digit:]]', names(lgcp.sim$raster), invert=TRUE, value=TRUE)
Sevents = grep("^events", names(lgcp.sim), value=TRUE)
#'


#+ aggregate

offsetCounty = tapply(
		kentuckyT$expected,
		list(kentuckyT$name2),
		sum)
kentucky$expected = offsetCounty[kentucky$County]
kentucky$logExpected = log(kentucky$expected)


for(Devent in Sevents) {
	projection(lgcp.sim[[Devent]]) = projection(kentucky)
	stuff = over(lgcp.sim[[Devent]], kentucky)$County
	lgcp.sim[[Devent]] = lgcp.sim[[Devent]][!is.na(stuff)] 
	stuff = table(factor(
					stuff, levels=unique(kentucky$County)
			))
	kentucky[[Devent]] = stuff[
			as.character(kentucky$County)
	]
}

kentuckyPoints = SpatialPoints(kentucky)
for(D in names(cov.ras)) {
	kentucky[[D]] = extract(
			cov.ras[[D]],		
			kentuckyPoints 	
	)
}
#'

#+ simPlot, fig.cap='truth', fig.subcap=c(Ssim, Sevents)

for(D in Ssim) {

	pCol = colourScale(
			lgcp.sim$raster[[D]],
			breaks=9, dec=-log10(0.25), col='Spectral',
			rev=TRUE, style='equal')
	map.new(lgcp.sim$raster,TRUE)
	plot(lgcp.sim$raster[[D]], 
			legend=FALSE,add=TRUE,
			col=pCol$col, breaks=pCol$breaks)
	plot(kMap,add=TRUE)
	legendBreaks("topleft", pCol)
}


for(D in Sevents) {
	map.new(lgcp.sim$raster,TRUE)
	plot(kMap,add=TRUE)
	plot(lgcp.sim[[D]],add=TRUE, col='#00FF0040')
}	

#'

#+ estimation, cache=TRUE, stuff=1
fit = list()
for(D in Sevents) {
	fit[[D]] = lgcp(
			data=lgcp.sim[[D]], 
      formula~w1+w2 + offset(offset),
			grid=squareRaster(kentucky, 80),
			border=kentucky,
      covariates=cov.ras,
      buffer=3,
      priorCI = list(sd=c(.1, 4),range=c(0.5,3)*100000)
	)
}
#'

#+ forResPLot
Splot = c('predict.exp', 'random.mean')
#'

#+ resPlot, fig.cap='posterior means', fig.subcap = rep(Splot, length(fit))

for(D in Splot) {
	for(Ds in names(fit)) {
	pCol = colourScale(
			fit[[Ds]]$raster[[D]],
			breaks=9, dec=1, col='Spectral',
			rev=TRUE, style='equal',
			opacity=0.7
			)
	map.new(fit[[1]]$raster, TRUE)
	plot(
			mask(fit[[Ds]]$raster[[D]],kentucky), 
			breaks=pCol$breaks, col=pCol$col, 
			add=TRUE, legend=FALSE)
	plot(kMap,add=TRUE)
	legendBreaks("right", pCol)
}
}
#'

#+ resTable, echo=FALSE
toPrint = lapply(fit, function(x)
			x$parameters$summary[,c(3,5)]
			)
names(toPrint) = gsub("events","s", names(toPrint))			
toPrint = do.call(cbind, toPrint)			
colnames(toPrint) = gsub("quant","q",colnames(toPrint))

knitr::kable(toPrint, digits=3)
#'



#+ estimationBym, cache=TRUE, stuff=4
fitBym = list()
for(D in Sevents) {
	kentucky$y = kentucky[[D]]
	fitBym[[D]] = bym(
			y ~ w1 + w2 + offset(logExpected), 
			data=kentucky, 
      priorCI = list(sd=c(.5, 0.1),propSpatial=c(0.5,0.1))
	)
}
#'



#+ resTableBym, echo=FALSE
toPrint = lapply(fit, function(x)
			x$parameters$summary[,c(3,5)]
)
names(toPrint) = gsub("events","s", names(toPrint))			
toPrint = do.call(cbind, toPrint)			
colnames(toPrint) = gsub("quant","q",colnames(toPrint))

knitr::kable(toPrint, digits=3)
#'



#+ ROC, cache=TRUE, stuff=2
res = spatialRoc(
		fit,
		rr=c(1,2,3), 
		truth = lgcp.sim,
		border=kentuckyT,
		random=FALSE
)

resR = spatialRoc(
		fit,
		rr=c(1,2,3), 
		truth = lgcp.sim,
		border=kentuckyT,
		random=TRUE
)

#'

#+ RocBym, cache=TRUE, stuff=4
resBym = spatialRoc(
		fitBym,
		rr=c(1,2,3), 
		truth = lgcp.sim,
		random=FALSE
)

resBymR = spatialRoc(
		fitBym,
		rr=c(1,2,3), 
		truth = lgcp.sim,
		border=NULL,
		random=TRUE
)

#'


#' # ROC



#+ ROCplot, fig.cap='ROC', fig.subcap=c('fitted lgcp','random lgcp', 'fitted bym', 'random bym'), fig.height=4, fig.width=4

		
Sspec = grep("spec", colnames(res), invert=TRUE, value=TRUE)
names(Sspec) = RColorBrewer::brewer.pal(length(Sspec), 'Dark2')
matplot(1-res[,'spec'], 
		res[,Sspec], 
		ylab='sens', xlab='1-spec',type='l',
		xlim=c(0,1), ylim=c(0,1), lty=1,
		col= names(Sspec)
)
legend("bottomright", col=names(Sspec), legend=Sspec, 
		lty=1, title='rr')




Sspec = grep("spec", colnames(resR), invert=TRUE, value=TRUE)
names(Sspec) = RColorBrewer::brewer.pal(length(Sspec), 'Dark2')
matplot(1-resR[,'spec'], 
		resR[,Sspec], 
		ylab='sens', xlab='1-spec',type='l',
		xlim=c(0,1), ylim=c(0,1), lty=1,
		col= names(Sspec)
)
legend("bottomright", col=names(Sspec), legend=Sspec, 
		lty=1, title='rr')

matplot(1-resBym[,'spec'], 
		resBym[,Sspec], 
		ylab='sens', xlab='1-spec',type='l',
		xlim=c(0,1), ylim=c(0,1), lty=1,
		col= names(Sspec)
)
legend("bottomright", col=names(Sspec), legend=Sspec, 
		lty=1, title='rr')

matplot(1-resBymR[,'spec'], 
		resBymR[,Sspec], 
		ylab='sens', xlab='1-spec',type='l',
		xlim=c(0,1), ylim=c(0,1), lty=1,
		col= names(Sspec)
)
legend("bottomright", col=names(Sspec), legend=Sspec, 
		lty=1, title='rr')

#'