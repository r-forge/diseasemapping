require('diseasemapping')
data('kentucky')
require('sp')

if(FALSE) {
	# must have an internet connection to do the following
	larynxRates= cancerRates("USA", year=1998:2002,site="Larynx")
	# get rid of under 10's
	larynxRates = larynxRates[-grep("_(0|5)$",names(larynxRates))]
	dput(larynxRates)
} else {
	larynxRates = structure(c(0, 0, 0, 0, 1e-06, 6e-06, 2.3e-05, 4.5e-05, 9.9e-05, 
					0.000163, 0.000243, 0.000299, 0.000343, 0.000308, 0.000291, 0.000217, 
					0, 0, 0, 1e-06, 1e-06, 3e-06, 8e-06, 1.3e-05, 2.3e-05, 3.5e-05, 
					5.8e-05, 6.8e-05, 7.5e-05, 5.5e-05, 4.1e-05, 3e-05), .Names = c("M_10", 
					"M_15", "M_20", "M_25", "M_30", "M_35", "M_40", "M_45", "M_50", 
					"M_55", "M_60", "M_65", "M_70", "M_75", "M_80", "M_85", "F_10", 
					"F_15", "F_20", "F_25", "F_30", "F_35", "F_40", "F_45", "F_50", 
					"F_55", "F_60", "F_65", "F_70", "F_75", "F_80", "F_85"))
	
}

kentucky = getSMR(kentucky, larynxRates, larynx,
		regionCode="County")

library('geostatsinla')
# this is in the examples

kBYM = bym(kentucky, observed ~ offset(logExpected) + poverty,
		priorCI = list(sdSpatial=c(0.1, 5), sdIndep=c(0.1, 5)))

# also try no covariate or prior

kBYM = bym(kentucky, observed ~ offset(logExpected))
 
if(!is.null(kBYM$parameters)) {
kBYM$par$summary

spplot(kBYM$data, "fitted.exp",
		at=seq(1.29, 1.36, len=6),
		col.regions= terrain.colors(5))

kBYM$data$exc1 = excProb(kBYM$inla$marginals.fitted.bym, log(1.2))


spplot(kBYM$data, "exc1",at=c(0, 0.2, 0.8,0.9, 1),
		col.regions=c('green','yellow','orange','red'))
}
# and try passing a data frame and adjacency matrix

library('spdep')
adjMat = poly2nb(kentucky, row.names =as.character(kentucky$County) )
kBYM = bym(kentucky@data, observed ~ offset(logExpected) + poverty,
		adjMat = adjMat, region.id="County",
		priorCI = list(sdSpatial=c(0.1, 5), sdIndep=c(0.1, 5)))
if(!is.null(kBYM$parameters)) {
kBYM$exc1 = excProb(kBYM$inla$marginals.fitted.bym, log(1.2))
}

# add subtract a few regions

kBYM = bym(kentucky@data[-(1:4),], observed ~ offset(logExpected) + poverty,
		adjMat = adjMat, region.id="County",
		priorCI = list(sdSpatial=c(0.1, 5), sdIndep=c(0.1, 5)))



if(!is.null(kBYM$parameters)) {
kBYM$exc1 = excProb(kBYM$inla$marginals.fitted.bym, log(1.2))
}

# intercept only, no offset


kBYM = bym(kentucky, observed ~ 1,
		priorCI = list(sdSpatial=c(0.1, 5), sdIndep=c(0.1, 5)))
if(!is.null(kBYM$parameters)) {
spplot(kBYM$data, "fitted.exp",
		at=c(0.8, 1.1, 1.2, 2, 4, 35),
		col.regions= terrain.colors(5))

kBYM$data$exc1 = excProb(kBYM$inla$marginals.fitted.bym, log(1.2))
}

# give spdf but some regions have no data
# but keep the 'county' column as is
kentucky@data[1:2,-grep("County", names(kentucky))] = NA 

kBYM = bym(kentucky, observed ~ offset(logExpected) + poverty,
		region.id="County",
		priorCI = list(sdSpatial=c(0.1, 5), sdIndep=c(0.1, 5)))
if(!is.null(kBYM$parameters)) {
kBYM$exc1 = excProb(kBYM$inla$marginals.fitted.bym, log(1.2))
}

