library("geostatsp")
data("swissRain")



# specify formula name of raster layer
swissFit = lgm(data=swissRain, formula=rain~ SRTM_1km,
		locations=80, covariates=swissAltitude,
		rough=1,  fixRough=TRUE, 
		boxcox=0.5, fixBoxcox=TRUE, 
		aniso=TRUE)	
names(swissFit)
swissFit$param


# specify formula using name of list element

swissFitAgain = lgm(data=swissRain, formula=rain~ elev,
		locations=80, covariates=list(elev=swissAltitude),
		rough=1,  fixRough=TRUE, 
		boxcox=0.5, fixBoxcox=TRUE, 
		aniso=TRUE)	
names(swissFitAgain)
swissFitAgain$param


swissFitAgain = lgm(data=swissRain, formula="rain",
		locations=80, covariates=swissAltitude,
		rough=1,  fixRough=TRUE, 
		boxcox=0.5, fixBoxcox=TRUE, 
		aniso=TRUE)	
names(swissFitAgain)
swissFitAgain$param
swissFitAgain = lgm(data=swissRain, formula="rain",
		locations=80, covariates=list(elev=swissAltitude),
		rough=1,  fixRough=TRUE, 
		boxcox=0.5, fixBoxcox=TRUE, 
		aniso=TRUE)	
names(swissFitAgain)
swissFitAgain$param

# land type, factor covariate
swissRes2 =  lgm(swissRain, locations=30, formula=rain ~ elev + factor(land),
		covariates=list(elev=swissAltitude,land=swissLandType), 
		boxcox=0.5, fixBoxcox=TRUE, 
		aniso=TRUE
)
swissRes2$summary


# simulated data (without a CRS)
# and all covariates are in 'data' object
myModel = c(intercept=0,variance=2^2,nugget=0.5^2, range=2.5,rough=2, 
		cov1=0.2, cov2=-0.5)
covariates = brick(
		xmn=0,ymn=0,xmx=10,ymx=10,
		ncols=200,nrows=200,nl=2)
values(covariates)[,1] = rep(seq(0,1,len=nrow(covariates)), ncol(covariates))
values(covariates)[,2] = rep(seq(0,1,len=nrow(covariates)), 
		rep(nrow(covariates), ncol(covariates)))
names(covariates) = c("cov1","cov2")

Npoints = 40
myPoints = SpatialPoints(cbind(runif(Npoints,0,10), runif(Npoints,0,10)))	
myPoints = SpatialPointsDataFrame(myPoints, 
		data=as.data.frame(extract(covariates, myPoints)))
myPoints$U = GaussRF(myPoints, param=myModel) 
myPoints$y= myModel["intercept"] +
		as.matrix(myPoints@data[,names(covariates)]) %*% 
		myModel[names(covariates)] +
		myPoints$U+
		rnorm(length(myPoints), 0, sqrt(myModel["nugget"]))

fitLikfit = likfitLgm(myPoints, trend=y~cov1+cov2, 
		param=c(range=1,nugget=0,rough=1)) 


Srange = seq(1, 3.5, len=10)
Slik = NULL
SlikWithN=NULL
Snugget=NULL
for(D in Srange) {
	Slik = c(Slik,
			loglikLgm(param=c(range=D,nugget=0,rough=1),
					data=myPoints, trend=y~cov1+cov2))
	temp = likfitLgm(paramToEstimate = "nugget",
			param=c(range=D,nugget=0,rough=1),
			data=myPoints, trend=y~cov1+cov2)
	SlikWithN = c(SlikWithN,
			temp$opt$value
			)	
	Snugget = c(Snugget, temp$opt$par["nugget"])		
}
plot(Srange, Slik)
plot(Srange, SlikWithN)
plot(Srange, Snugget)



# run lgm without providing covariates
fitMLE =  lgm(myPoints, locations=10, formula=y~ cov1 + cov2, 
		rough=1, fixRough=TRUE)

c(fitMLE$summary["range","estimate"], fitLikfit$summary["range","estimate"])


# now give covariates as raster brick
fitMLE =  lgm(myPoints, locations=10, formula=y~ cov1 + cov2, 
		covariates=covariates,
		rough=1, fixRough=TRUE)
c(fitMLE$summary["range","estimate"], fitLikfit$summary["range","estimate"])

# now give covariates as list
fitMLE =  lgm(myPoints, locations=10, formula=y~ cov1 + cov2, 
		covariates=list(cov1=covariates[["cov1"]],
				cov2 = covariates[["cov2"]]),
		rough=1, fixRough=TRUE)
c(fitMLE$summary["range","estimate"], fitLikfit$summary["range","estimate"])


# not remove covariates from data
myPoints = SpatialPointsDataFrame(SpatialPoints(myPoints),
		data=myPoints@data[,"y",drop=FALSE])

# now give covariates as raster brick
fitMLE =  lgm(myPoints, locations=10, formula=y~ cov1 + cov2, 
		covariates=covariates,
		rough=1, fixRough=TRUE)
c(fitMLE$summary["range","estimate"], fitLikfit$summary["range","estimate"])


