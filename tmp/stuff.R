
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