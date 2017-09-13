# TO DO: fix x > 2
# fix large matrices

options(width=200)

Rcpp::compileAttributes("/home/patrick/workspace/diseasemapping/pkg/geostatsgpu")
devtools::install("/home/patrick/workspace/diseasemapping/pkg/geostatsgpu")
#devtools::load_all("/home/patrick/workspace/diseasemapping/pkg/geostatsgpu")
#devtools::load_all("/home/patrick/workspace/diseasemapping/pkg/geostatsp")
devtools::reload("/home/patrick/workspace/diseasemapping/pkg/geostatsgpu")
library('geostatsgpu')
library('geostatsp')

x = as.matrix(expand.grid(seq(1,by=1,len=23), seq(201,by=1,len=14)))
coordsV = gpuR::vclMatrix(x)
D3 <- vclMatrix(
    data=-1, 
    nrow(coordsV), nrow(coordsV),
    type='double'
)

myParams = c(shape=4.5, range=1.5, variance = 2, nugget = 0, anisoRatio = 2, anisoAngleRadians = pi/4)
myParams2 = c(shape=4.5, range=1.5, variance = 2, nugget = 0, anisoRatio = 1, anisoAngleRadians = -pi/4)

system.time(v1 <- distGpu(coordsV, myParams, D3, 'variance'))

mmat=as.matrix(geostatsp::matern(sp::SpatialPoints(as.matrix(coordsV)), myParams))
xxmat = log(sqrt(8*myParams['shape'])*as.matrix(dist(x))/myParams['range'])
signif(
    data.frame(
a=      as.matrix(v1)[,1],
b=      xxmat[,1],
c=      mmat[,1]
)[1:17,], 2)

plot(xxmat[,1], log(abs(as.matrix(v1)[,1]/mmat[,1])),
    ylim = c(-1,1)/12)


round(as.matrix(v1) - as.matrix(dist(x))^2, digits=5)[1:12,1:12]

as.matrix(dist(x))[1:5,1:24]

#myDist = outer(x[,1] + 1i*x[,2],x[,1] + 1i*x[,2],FUN='-')
#Mod(myDist)^2

library('geostatsp')
system.time(v2 <- geostatsp::matern(
    sp::SpatialPoints(as.matrix(coordsV)),
    myParams
))
as.matrix(v2)[1:11,1:4]


cbind(as.matrix(v1)[,1], v2[,1])

as.matrix(v1) - v2
as.matrix(v1) - Mod(myDist)^2


options(width=120)
v1[1:7,1:7]
v2[1:7,1:7]


(Mod(myDist)^2)[11,1]


(v2 - as.matrix(v1))[1:7,1:7]
(v2 - as.matrix(v1))[seq(by=1, len=7, to=dim(v2)[1]), seq(by=1, len=7, to=dim(v2)[2])]


v2[1:5,1:4]
as.matrix(v1)[1:5,1:4]


quantile((as.matrix(v1) - v2)[lower.tri(v2)])



# development version of gpuR

cl_args <- setup_opencl(objects = c(
        "scalar",
        "vclVector", 
        "vclVector"),
    intents = c("IN", "OUT", "IN"),
    queues = list("SAXPY", "SAXPY", "SAXPY"),
    kernel_maps = c("x", "y", "a"))



