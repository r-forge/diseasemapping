library('mapmisc')
Spackages = c('rgdal', 'rgeos', 'geosphere', 'maptools')

if(all(unlist(mapply(requireNamespace, package=Spackages, MoreArgs=list(quietly=TRUE))))){
	library('rgdal')

	data('wrld_simpl', package='maptools')
	
	country='Japan'
	Dcountry  = grep(country, wrld_simpl$NAME)
	x=wrld_simpl[Dcountry,]

	myCrsO = ocea(x, angle=25)
	xTcrop = wrapPoly(wrld_simpl, myCrsO)
	DcountryT  = grep(country, xTcrop$NAME)
	
	map.new(xTcrop)
	plot(xTcrop,add=TRUE)
	plot(xTcrop[DcountryT,], col='red', add=TRUE)
	rgdal::llgridlines(xTcrop, col='blue', 
			easts=seq(-180,180,by=30),
			norths=seq(-90,90,by=30))
	plot(attributes(myCrsO)$circleTrans, add=TRUE, col='#00FF0020', pch=1)
	
	country='Madagascar'
	Dcountry  = grep(country, wrld_simpl$NAME)
	x=wrld_simpl[Dcountry,]
	
	myCrsMoll = moll(x, oblique=TRUE)
	xTcrop = wrapPoly(wrld_simpl, myCrsMoll)
	DcountryT  = grep(country, xTcrop$NAME)
	
	map.new(xTcrop)
	plot(xTcrop,add=TRUE)
	plot(xTcrop[DcountryT,], col='red', add=TRUE)
	rgdal::llgridlines(xTcrop, col='blue', 
			easts=seq(-180,180,by=30),
			norths=0#seq(-90,90,by=30)
	)
	
}