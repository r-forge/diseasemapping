
x = raster(nrows=11,ncols=11,xmn=0,xmx=10,ymn=0,ymx=10)

temp=GaussRF(x, model="whittle", 
		param=c(mean=0, variance=1, nugget=0, 
				scale=2, alpha=2))


x = SpatialPoints(cbind(1:4, 11:14))

 temp=GaussRF(x, model="whittle", 
		param=c(mean=0, variance=1, nugget=0, 
				scale=2, alpha=2))


x = cbind(1:4, 1:4)
temp=GaussRF(x, model="whittle", 
		param=c(mean=0, variance=1, nugget=0, 
				scale=2, alpha=2))

Ncell = 25

# as in example
require('geostatsp')

data('swissRain')
swissRain$lograin = log(swissRain$rain)
debug(glgm)
	swissFit =  glgm(lograin ~ CHE_alt, swissRain, Ncell, 
			#covariates=swissAltitude, 
			family="gaussian", buffer=20000,
			priorCI=list(sd=c(0.2, 2), range=c(50000,500000)), 
			control.mode=list(theta=c(1.9,0.15,2.6),restart=TRUE),
			control.family=list(hyper=list(prec=list(prior="loggamma", param=c(.1, .1))))
	)
	
	
	swissFit = lgm(data=swissRain, formula=rain~ CHE_alt,
			grid=80, #covariates=swissAltitude,
			shape=1,  fixShape=TRUE, 
			boxcox=0.5, fixBoxcox=TRUE, 
			aniso=TRUE)	
	
	
	list(mal_data)
	
	names(mal)
	
	library(mapmisc)
	malmap = openmap(mal)

	map.new(mal)
	plot(malmap,add=TRUE)
	plot(mal,add=TRUE)
	
	
#altitude
	map.new(mal)
	plot(malmap,add=TRUE)
	plot(ALT, col=terrain.colors(100), main="elevation", alpha=0.5, add=TRUE)
	points(mal)
	
	
#NDVI
map.new(mal)
plot(malmap,add=TRUE)
	plot(NDVI, main="mean ndvi",alpha=0.5,col=terrain.colors(100), add=TRUE)
	points(mal)
	
#Population
map.new(mal)
plot(malmap,add=TRUE)
plot(POP, main="Population Density",alpha=0.5,col=terrain.colors(100), add=TRUE)
	points(mal)
	
	
	
#Prepare data for modeling
	covList = list(alt=ALT, maxtempwmq=MAXTWMM, precwetmn=PRECWTM,
			precdrm=PRECDRM, #precwetq=PRECWQ, 
			preccoldq=PRECCQ, 
			ndvi=NDVI, pop=POP)
	
for(D in names(covList)) {
	map.new(mal)
	plot(malmap,add=TRUE)
	plot(POP, alpha=0.5,col=terrain.colors(100), add=TRUE)
	points(mal, col='#FF000070')
	mtext(D,side=3, adj=0.8,line=-3,cex=2)
	invisible(readline(prompt="Press [enter] to continue"))	
}	





#Run the model
	malFit<-glgm(
			formula = POS ~ alt + maxtempwmq + precwetmn + 
					precdrm + #precwetq + 
					preccoldq + ndvi + pop, 
			data = mal, grid = 50,
			covariates = covList, family = "binomial", 
			Ntrials = mal$TOT, shape = 1,buffer = 50000, 
			priorCI = list(sd = c(0.2, 4), range = c(20000, 5e+05)))
	
	grid=50;data=mal;buffer=50000;covariates=covList
	grid=gridRaster = squareRaster(data, grid)
	cellsBoth = cellsBuffer(grid, buffer)			
	cellsSmall = cellsBoth$small
	
	rmethod = rep("bilinear", length(names(covariates)))
	names(rmethod) = names(covariates)
	rmethod['precwetq'] = "ngb"
	covariatesStack = stackRasterList(covariates, 
			template=cellsSmall, 
			method=rmethod)
	