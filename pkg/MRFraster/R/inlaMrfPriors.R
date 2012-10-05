inlaMrfPriors <-
function(rangePmode=1, rangeShape=1, sdPmode=0.1, varShape=2, 
		gridsize=1, cellVar="cellID") {
	
	if( length(grep('^Raster',class(gridsize)))  ) {
		gridsize = xres(gridsize)
	}
	
	# specify prior mode of range in km, range is inverse gamma distribution
	# INLA wants prior on gridsize/range to be gamma
	range = list(prior='loggamma',
			param = c(shape=rangeShape, scale=(rangeShape+1)*rangePmode/gridsize)
	)
#	range$initial=log(range$param["shape"]/range$param["scale"])
	
	attributes(range)$ci95RangeKm = sort(gridsize/qgamma(c(0.025, 0.975), 
					shape=range$param["shape"], scale=1/range$param["scale"]))
	
	prec = list(prior='loggamma',
			param=c(shape=varShape, scale=sdPmode^2*(varShape+1))
	)
#	prec$initial=log(prec$param["shape"]/prec$param["scale"])
	
	attributes(prec)$ci95sd = sort(1/sqrt(qgamma(c(0.025, 0.975), 
							shape=prec$param["shape"], scale=1/prec$param["scale"])))
	
	return(list(range=range, prec=prec, 
				formula=
					as.formula(paste(".~.+f(",cellVar,", model='matern2d',nrow=",
									attributes(dataForInla)$nrow, 
									", ncol=",
									attributes(dataForInla)$ncol, 
									", nu=2,hyper = list(range=",
									as.character(thePriors['range']),
									", prec=",
									as.character(thePriors['prec']),
									") )")) 
	))
}
