library("geostatsp")

param = c(range=0, rough=1.5,	aniso.ratio=2, aniso.angle.degrees=-25)

matern(c(0, 0.001, 100000), param=param)

# example with raster
myraster = raster(nrows=40,ncols=60,xmn=-3,xmx=3,ymn=-2,ymx=2)

# plot correlation of each cell with the origin
myMatern = matern(myraster, c(0,0), param=param)
myMatern[1:3,1:3]

plot(myMatern)


# correlation matrix for all cells with each other
myraster = raster(nrows=4,ncols=6,xmn=-3,xmx=3,ymn=-2,ymx=2)
myMatern = matern(myraster, param=c(range=0, rough=2))
dim(myMatern)
myMatern[1:3,1:3]

mypoints = SpatialPointsDataFrame(cbind(runif(5), runif(5)),data=data.frame(id=1:5))
matern(mypoints, param=param)
