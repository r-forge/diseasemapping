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
output.ras = squareRaster(kentuckyT, 300)
cov.ras = list()
cov.ras$w1 = cov.ras$w2 = output.ras
values(cov.ras$w1) = yFromCell(output.ras, 1:ncell(output.ras))/500000
values(cov.ras$w2) = xFromCell(output.ras, 1:ncell(output.ras))/500000
#'


#+ offsets, cache=TRUE, stuff=1
larynxRates= cancerRates("USA", year=1998:2002,site="Larynx")
larynxRates = larynxRates * 5
kentuckyT = getSMR(kentuckyT, larynxRates, regionCode="County")

kentuckyOffset = spdfToBrick(
		x=kentuckyT[,'expected'],
		template=output.ras
		)
cov.ras$kentuckyOffset = kentuckyOffset
cov.ras$offset = log(kentuckyOffset) 
#'

#+ simB, cache=TRUE, stuff=1
U = c(mean = 0, 
		variance=0.25^2, 
		range=80000,
		shape=2)

set.seed(0)
lgcp.sim = simLgcp(
		param=U,
		covariates=cov.ras,
		betas=c(w1=-0.1,w2=0.025),
		rasterTemplate=output.ras, 
		n=12,
		offset='offset')
#'

#+ forSimPlot
Ssim = grep('^linear|^intensity|^offset|^w[[:digit:]]', 
		names(lgcp.sim$raster), invert=TRUE, value=TRUE)
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


#+ eventsPlot, fig.cap='events',  fig.subcap=Sevents

for(D in Sevents) {
	map.new(kentucky,TRUE)
	plot(kMap,add=TRUE)
	plot(lgcp.sim[[D]],add=TRUE, col='#00FF0040')
}	

#'

#+ estimationGeostatsp
geostatspFile = 'geostatsp.RData'
if(!file.exists(geostatspFile)) {
fitLgcp = list()
for(D in Sevents) {
	fitLgcp[[D]] = lgcp(
			data=lgcp.sim[[D]], 
      formula=~w1+w2 + offset(kentuckyOffset, log=TRUE),
			grid=squareRaster(kentucky, 100),
      covariates=cov.ras,
      buffer=3,
      priorCI = list(sd=c(u=0.25, alpha=0.1),range=c(0.5,3)*100000)
	)
}
save(fitLgcp, file=geostatspFile)
} else {
	load(geostatspFile)
}
#'


#+ estimationBym
bymFile = 'bym.RData'
if(!file.exists(bymFile)){
	
fitBym = list()
for(D in Sevents) {
	kentucky$y = kentucky[[D]]
	fitBym[[D]] = bym(
			y ~ w1 + w2 + offset(logExpected), 
			data=kentucky, 
      priorCI = list(sd=c(.5, 0.1),propSpatial=c(0.5,0.1))
	)
}
save(fitBym, file=bymFile)
} else {
	load(bymFile)
}
#'




#+ forResPLot
Splot = c('predict.exp', 'random.mean')
SplotMat = expand.grid(var=Splot, sim=names(fitLgcp))
SplotSubcap = as.vector(paste(SplotMat[,'var'], gsub("events", "", SplotMat[,'sim'])))
SplotBym = c('fitted.exp', 'random.mean')
SplotMatBym = expand.grid(var=SplotBym, sim=names(fitLgcp))
#'

#+ resPlot, fig.cap='posterior means lgcp', fig.subcap = SplotSubcap
for(Drow in 1:nrow(SplotMat)) {

	D = as.character(SplotMat[Drow,'var'])
	Ds = as.character(SplotMat[Drow,'sim'])
	pCol = colourScale(
			fitLgcp[[Ds]]$raster[[D]],
			breaks=9, dec=1, col='Spectral',
			rev=TRUE, style='equal'
			)
	map.new(kentucky, TRUE)
	plot(
			mask(fitLgcp[[Ds]]$raster[[D]],kentucky), 
			breaks=pCol$breaks, col=pCol$col, 
			add=TRUE, legend=FALSE)
	plot(kMap,add=TRUE)
	legendBreaks("right", pCol)
}
#'

#+ resTable, echo=FALSE
toPrint = lapply(fitLgcp, function(x)
			x$parameters$summary[,c(3,5)]
			)
names(toPrint) = gsub("events","s", names(toPrint))			
toPrint = do.call(cbind, toPrint)			
colnames(toPrint) = gsub("quant","q",colnames(toPrint))

knitr::kable(toPrint, digits=3)
#'





#+ resPlotB, fig.cap='posterior means BYM', fig.subcap = SplotSubcap

for(Drow in 1:nrow(SplotMatBym)) {
	
	D = as.character(SplotMatBym[Drow,'var'])
	Ds = as.character(SplotMatBym[Drow,'sim'])

	pCol = colourScale(
				fitBym[[Ds]]$data[[D]],
				breaks=9, dec=2, col='Spectral',
				rev=TRUE, style='equal'
		)
		map.new(kentucky, TRUE)
		plot(
				fitBym[[Ds]]$data, 
				col=pCol$plot, 
				add=TRUE, border='#00000030')
		plot(kMap,add=TRUE)
		legendBreaks("right", pCol)
}
#'

#+ resTableBym, echo=FALSE
toPrint = lapply(fitBym, function(x)
			x$parameters$summary[,c(3,5)]
)
names(toPrint) = gsub("events","s", names(toPrint))			
toPrint = do.call(cbind, toPrint)			
colnames(toPrint) = gsub("quant","q",colnames(toPrint))

knitr::kable(toPrint, digits=3)
#'



#+ ROC, cache=TRUE, stuff=1

theRR = c(1.05, 1.2, 1.5)

resLgcp = spatialRoc(
		fitLgcp,
		rr=theRR, 
		truth = lgcp.sim,
		border=kentuckyT,
		random=FALSE
)

resLgcpR = spatialRoc(
		fitLgcp,
		rr=theRR, 
		truth = lgcp.sim,
		border=kentuckyT,
		random=TRUE
)


resBym = spatialRoc(
		fit=fitBym,
		rr=theRR, 
		prob = 2^seq(-1, -10, len=10),
		truth = lgcp.sim,
		random=FALSE
)

resBymR = spatialRoc(
		fitBym,
		rr=theRR, 
		truth = lgcp.sim,
		random=TRUE
)

#'


#' # ROC

#+ junk, eval=FALSE, include=FALSE

sum(values(kentuckyOffset),na.rm=TRUE)*prod(res(kentuckyOffset))
sum(values(cov.ras$offset),na.rm=TRUE)*prod(res(cov.ras$offset))
sum(values(covariates$offset),na.rm=TRUE)*prod(res(covariates$offset))
sum(kentuckyT$expected)
sum(fitLgcp[[1]]$inla$.args$data$offset*prod(res(fitLgcp[[1]]$raster)),na.rm=TRUE)
sum(exp(fitLgcp[[1]]$inla$.args$data$logoffset+
						fitLgcp[[1]]$inla$.args$data$logCellSize))

fit=fitLgcp

fit=fitBym
rr=theRR
truth = lgcp.sim
border=NULL
random=FALSE 
prob = NULL
spec = seq(0,1,by=0.01)

Dbin = '3'

Dsim = names(marginals)[1]


stuff = values(truthCut[[3]])[which(values(templateID==4))]
stuff
truthCdf[64,,3]
Dlevel = '3.5'
truthOver[64,,3]
truthUnder[64,,3]
pMat[64,]
predOver[64,]
tpMat[64,]
fpMat[64,]

Dbin
Dmidpoint
Dsim

predOver[64,]*truthOver[64,Dmidpoint,Dsim]
(1-predOver[64,])*truthUnder[64,Dmidpoint,Dsim]

tpMat[64,]
fpMat[64,]

res$tP[,3,3]
res$tN[,3,3]


plot(truthCut[[3]])
plot(kentucky[4,],add=TRUE, border='red')

ob = attributes(resBym)$orig
sb = attributes(resBym)$sim
sl = attributes(resLgcp)$sim
ol = attributes(resLgcp)$orig

Dbin = '3'
Dsim = 1
plot(sb[,Dbin,Dsim,],  type='o', xlim=c(0,0.1), ylim=c(0,0.8))
lines(sl[,Dbin,Dsim,], col='red', type='o')
#text(sl[,Dbin,3,], substr(dimnames(sl)[[1]],1,6), cex=0.9, col='blue')
text(sb[,Dbin,Dsim,], substr(dimnames(sb)[[1]],1,4), cex=0.8)

sb[,3,3,]


sb[seq(to=max(which(sb[,3,3,'onemspec']>0))+1,by=1,len=6),3,3,]
sl[seq(to=max(which(sl[,3,3,'onemspec']>0))+1,by=1,len=6),3,3,]

ob['0.5',,]
ol['0.5',,]

sl['0.5',,1,]
sb[15,,1,]

#'


#+ ROCplot, fig.cap='ROC', fig.subcap=c('fitted','random'), fig.height=5, fig.width=6, fig.ncol=1

		
Sspec = grep("spec", colnames(resLgcp), invert=TRUE, value=TRUE)
names(Sspec) = RColorBrewer::brewer.pal(length(Sspec), 'Dark2')

matplot(resLgcp[,'onemspec'], 
		resLgcp[,Sspec], 
		ylab='sens', xlab='1-spec',type='l',
		xlim=c(0,0.5), ylim=c(0,1), lty=1,
		col= names(Sspec)
)

matlines(resBym[,'onemspec'], 
		resBym[,Sspec], 
		lty=3, lwd=4,
		col= names(Sspec)
)

legend("bottomright", col=names(Sspec), legend=Sspec, 
		lty=1, title='rr')
legend("right", lty=c(1,3), legend=c('lgcp','bym'))


matplot(resLgcpR[,'onemspec'], 
		resLgcpR[,Sspec], 
		ylab='sens', xlab='1-spec',type='l',
		xlim=c(0,0.5), ylim=c(0,1), lty=1,
		col= names(Sspec)
)


matlines(resBymR[,'onemspec'], 
		resBymR[,Sspec], 
		lty=3, lwd=4,
		col= names(Sspec)
)
legend("bottomright", col=names(Sspec), legend=Sspec, 
		lty=1, title='rr')
legend("right", lty=c(1,3), legend=c('lgcp','bym'))

#'