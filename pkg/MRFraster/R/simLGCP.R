simLGCP <-
function(N=101,mu=3, beta=c(1,1.5), sigma=1, 
		range=0.3, roughness=2){
library(RandomFields)
library(raster)

model="whittle"

# try mu 3 and 4
param = list(mu=mu, beta=beta, sigma=sigma, 
		range=range, roughness=roughness,model=model)

mybbox = cbind(c(0,0), c(1,1))
dimnames(mybbox)=list(c("x","y"), c("min","max"))

#myextent = extent(mybbox["x","min"], mybbox["x","max"],mybbox["y","min"], mybbox["y","max"])
#rasterTemplate = raster(nrows=N, ncols=N, ext=myextent)



x <- seq(mybbox["x","min"], mybbox["x","max"], len=N) 
y <- seq(mybbox["y","min"], mybbox["y","max"], len=N) 

U <- GaussRF(x=x, y=y, model=param$model, grid=TRUE,
        param=c(mean=0, variance=1, nugget=0, 
                scale=param$range, nu=param$roughness))

U=raster(matrixForRaster(U),
		xmn=mybbox["x","min"], xmx=mybbox["x","max"],
		ymn=mybbox["y","min"], ymx=mybbox["y","max"])
U@layernames='U'


U = param$sigma*(U-mean(U[]))/sd(U[]) 


offset = outer(x, 1i*y, FUN='+')
offset = offset - mean(x)-1i*mean(y)
offset = 2 - Mod(offset)

offset = matrix(offset, 
		length(x), length(y)) 
offset = raster(matrixForRaster(offset), 
	xmn=mybbox["x","min"], xmx=mybbox["x","max"],
	ymn=mybbox["y","min"], ymx=mybbox["y","max"])
offset@layernames='offset'


covs = stack(raster(matrix(x, length(x), length(y)), 
		xmn=mybbox["x","min"], xmx=mybbox["x","max"],
		ymn=mybbox["y","min"], ymx=mybbox["y","max"])
, 
raster(matrix(y, length(x), length(y), byrow=T), 
		xmn=mybbox["x","min"], xmx=mybbox["x","max"],
		ymn=mybbox["y","min"], ymx=mybbox["y","max"])
)
covs@layernames=c('x','y')


lambda = offset + U + param$mu + param$beta[1]*covs[[1]]+ 
		param$beta[2]*covs[[2]]
lambda = exp(lambda)
lambda@layernames = 'lambda'

pointsRast = lambda
pointsRast[] = rpois(length(lambda), xres(lambda)*yres(lambda)*lambda[])

toRep = rep(1:length(pointsRast), pointsRast[])

pointsSP = SpatialPoints(pointsRast)[toRep]

pointsSP@bbox = mybbox

# jitter points within grid cells
	xres2 = xres(lambda)/2
	yres2 = yres(lambda)/2
	pointsSP@coords = pointsSP@coords + 
			cbind(runif(length(pointsSP), -xres2, xres2),
					runif(length(pointsSP), -xres2, xres2)
					)



return(list(points=pointsSP, params=param, 
				lambda=lambda/exp(offset),
				intensity=lambda,
				U=U, covs = covs,
				offset=offset
) )






}
