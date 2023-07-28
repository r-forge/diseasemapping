havePackages = c(
  'INLA' = requireNamespace('INLA', quietly=TRUE)
)

print(havePackages)

library('diseasemapping')
data('kentucky')
kentucky = unwrap(kentucky)





kentucky = getSMR(kentucky, larynxRates, larynx,
		regionCode="County")

if(all(havePackages)){

  kBYM = bym(
			formula = observed ~ offset(logExpected) + poverty,
      data=values(kentucky),
			adjMat = terra::adjacent(kentucky),
      prior = list(sd=c(0.1, 5), propSpatial=c(0.1, 5)),
			region.id="County",
			control.predictor=list(compute=TRUE)
  )
	kBYM$parameters$summary

	pdf("priorPostKentucky.pdf")
	plot(kBYM$parameters$sdSpatial$posterior, type='l', 
			xlim=c(0,max(kBYM$parameters$sdSpatial$priorCI)))
	lines(kBYM$parameters$sdSpatial$prior, col='blue')
	legend('topright', lty=1, col=c('black','blue'), legend=c('posterior','prior'))
	dev.off()

	
	


# also try no covariate or prior

kBYM = bym(
		formula = observed ~ offset(logExpected),
		data=kentucky)


if(require('geostatsp', quietly=TRUE)) {
 	kBYM$data$exc1 = geostatsp::excProb(kBYM$inla$marginals.fitted.bym, log(1.2))
} else {
	kBYM$data$exc1 = rep(NA, length(kBYM$data))
}

kBYM$par$summary

if(require('mapmisc', quietly=TRUE)) {

colFit = colourScale(kBYM$data$fitted.exp,
		breaks=6, dec=3)
	
plot(kBYM$data, col=colFit$plot)
legendBreaks('topleft', colFit)

colExc = colourScale(kBYM$data$exc1 ,
		style='fixed',
		breaks=c(0, 0.2, 0.8,0.9, 1), 
		col=rainbow, rev=TRUE
	)

	plot(kBYM$data, col=colExc$plot)
	legendBreaks('topleft', colExc)
 		
}
# and try passing a data frame and adjacency matrix

	
adjMat = spdep::poly2nb(kentucky, row.names =as.character(kentucky$County) )
kBYM = bym(data=kentucky@data, formula=observed ~ offset(logExpected) + poverty,
		adjMat = adjMat, region.id="County",
		priorCI = list(sdSpatial=c(0.1, 5), sdIndep=c(0.1, 5)))

kBYM$par$summary

# subtract a few regions

kBYM = bym(
    formula=observed ~ offset(logExpected) + poverty,
    data=kentucky@data[-(1:4),],  
 	  adjMat = adjMat, region.id="County",
		priorCI = list(sdSpatial=c(0.1, 5), sdIndep=c(0.1, 5)))
 

kBYM$par$summary

# intercept only, no offset


kBYM = bym(data=kentucky,  formula=observed ~ 1,
		priorCI = list(sdSpatial=c(0.1, 5), sdIndep=c(0.1, 5)))

kBYM$par$summary


if(require('mapmisc', quietly=TRUE)) {
	
	colFit = colourScale(kBYM$data$fitted.exp,
			breaks=6, dec=1)
	
	plot(kBYM$data, col=colFit$plot)
	legendBreaks('topleft', colFit)
	
}

 

# give spdf but some regions have no data
# but keep the 'county' column as is
kentucky@data[1:2,-grep("County", names(kentucky))] = NA 

kBYM = bym(observed ~ offset(logExpected) + poverty,
		kentucky, 
		region.id="County",
		priorCI = list(sdSpatial=c(0.1, 5), sdIndep=c(0.1, 5)))

 
kBYM$par$summary


# missing value in a categorical variable

pCuts = quantile(kentucky$poverty, na.rm=TRUE)
kentucky$povertyFac = cut(kentucky$poverty, 
		breaks = pCuts,
		labels = letters[seq(1,length(pCuts)-1)])
kentucky$povertyFac[c(2,34,100)] = NA

kBYM = bym(
		formula = observed ~ offset(logExpected) + povertyFac,
		data = kentucky, 
		region.id="County",
		priorCI = list(sdSpatial=c(0.1, 5), sdIndep=c(0.1, 5))
)


kBYM$par$summary
}


