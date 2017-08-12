
devtools::load_all("/home/patrick/workspace/diseasemapping/pkg/geostatsgpu")

library(gpuR)
x = as.matrix(expand.grid(seq(1,by=1,len=41), seq(1,by=1,len=41)))
coordsV = gpuR::vclMatrix(x)
D3 <- vclMatrix(
    data=-0.1, 
    nrow(coordsV), nrow(coordsV),
    type='double'
)

myParams = c(shape=1.5, range=1, variance = 1, nugget = 0, anisoRatio = 1, anisoAngleRadians = pi/4)

v1 = maternGpu(coordsV, myParams, D3, 'variance')

as.matrix(v1)[1,]

as.matrix(v1)[,1]

as.matrix(v1)[1:12,1:8]

v2 = geostatsp::matern(
    sp::SpatialPoints(as.matrix(coordsV)),
    myParams
)

options(width=120)
(v2 - as.matrix(v1))[1:7,1:7]
(v2 - as.matrix(v1))[seq(by=1, len=7, to=dim(v2)[1]), seq(by=1, len=7, to=dim(v2)[2])]


v2[1:5,1:4]
as.matrix(v1)[1:5,1:4]


quantile((as.matrix(v1) - v2)[lower.tri(v2)])

