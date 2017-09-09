
# TO DO: fix x > 2
# fix large matrices

devtools::load_all("/home/patrick/workspace/diseasemapping/pkg/geostatsgpu")

library(gpuR)
x = as.matrix(expand.grid(seq(1,by=1,len=4), seq(201,by=1,len=5)))
coordsV = gpuR::vclMatrix(x)
D3 <- vclMatrix(
    data=-0.1, 
    nrow(coordsV), nrow(coordsV),
    type='double'
)

myParams = c(shape=1.5, range=2, variance = 1, nugget = 0, anisoRatio = 1, anisoAngleRadians = pi/4)

system.time(v1 <- maternGpu(coordsV, myParams, D3, 'variance'))
options(width=200)
as.matrix(v1)[1:11,1:4]

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