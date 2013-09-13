
# as in example
require('geostatsinla')
require('sp')
data('swissRain')
swissRain$lograin = log(swissRain$rain)
swissFit =  glgm(swissRain, cells=30, formula="lograin",
		covariates=swissAltitude, family="gaussian", buffer=20000,
		priorCI=list(sd=c(0.2, 2), range=c(50000,500000)), 
		control.mode=list(theta=c(1.9,0.15,2.6),restart=TRUE),
		control.family=list(hyper=list(prec=list(prior="loggamma", param=c(.1, .1))))
)

swissFit$parameters$summary

swissExc = excProb(swissFit$inla$marginals.random$space, 0, swissFit$raster)
plot(swissExc, breaks = c(0, 0.2, 0.8, 0.95, 1.00001), 
		col=c('green','yellow','orange','red'))	
plot(swissBorder, add=TRUE)		


# intercept only
swissFit =  glgm(swissRain, cells=30, formula=lograin~1,
		covariates=swissAltitude, family="gaussian", buffer=20000,
		priorCI=list(sd=c(0.2, 2), range=c(50000,500000)), 
		control.mode=list(theta=c(1.9,0.15,2.6),restart=TRUE),
		control.family=list(hyper=list(prec=list(prior="loggamma", param=c(.1, .1))))
)

swissFit$parameters$summary

swissExc = excProb(swissFit$inla$marginals.random$space, 0, swissFit$raster)
plot(swissExc, breaks = c(0, 0.2, 0.8, 0.95, 1.00001), 
		col=c('green','yellow','orange','red'))	
plot(swissBorder, add=TRUE)		


# now with formula
swissFit =  glgm(swissRain, cells=30, 
		formula=lograin~ SRTM_1km,
		covariates=swissAltitude, family="gaussian", buffer=20000,
		priorCI=list(sd=c(0.2, 2), range=c(50000,500000)), 
		control.mode=list(theta=c(1.9,0.15,2.6),restart=TRUE),
		control.family=list(hyper=list(prec=list(prior="loggamma", 
								param=c(.1, .1))))
)

# formula, named list elements
swissFit =  glgm(swissRain, cells=30, 
		formula=lograin~ elev,
		covariates=list(elev=swissAltitude), 
		family="gaussian", buffer=20000,
		priorCI=list(sd=c(0.2, 2), range=c(50000,500000)), 
		control.mode=list(theta=c(1.9,0.15,2.6),restart=TRUE),
		control.family=list(hyper=list(prec=list(prior="loggamma", 
								param=c(.1, .1))))
)

# categorical covariates
swissFit =  glgm(swissRain, cells=30, 
formula=rain ~ elev + factor(land),
covariates=list(elev=swissAltitude,land=swissLandType), 
family="gaussian", buffer=20000,
priorCI=list(sd=c(0.2, 2), range=c(50000,500000)), 
control.mode=list(theta=c(1.9,0.15,2.6),restart=TRUE),
control.family=list(hyper=list(prec=list(prior="loggamma", 
						param=c(.1, .1))))
)


# put some missing values in covaritates
temp = values(swissAltitude)
temp[seq(10000,12000)] = NA
values(swissAltitude) = temp
swissFit =  glgm(swissRain, cells=30, 
		formula=rain ~ elev + factor(land),
		covariates=list(elev=swissAltitude,land=swissLandType), 
		family="gaussian", buffer=20000,
		priorCI=list(sd=c(0.2, 2), range=c(50000,500000)), 
		control.mode=list(theta=c(1.9,0.15,2.6),restart=TRUE),
		control.family=list(hyper=list(prec=list(prior="loggamma", 
								param=c(.1, .1))))
)




