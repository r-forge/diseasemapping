library(gpuR)

coords = as.matrix(expand.grid(1:5, 1:5))/5
storage.mode(coords) = 'numeric'
coordsV = vclMatrix(coords)

distV = dist(coordsV, 'sqEuclidean')

x=coords
D1 <- vclMatrix(as.double(0), nrow=nrow(x), ncol=nrow(x), type='double', ctx_id=coordsV@.context_index)

gpuR:::vclMatrix_euclidean(
    A=coordsV, 
    D=D1,
    diag = TRUE,
    upper = FALSE,
    p = 2,
    squareDist = TRUE)

D2 <- vclMatrix(as.double(0),
    nrow=nrow(x), ncol=nrow(x), type='double', 
    ctx_id=coordsV@.context_index)

gpuR:::cpp_vclMatrix_eucl(
    coordsV@address, 
    D2@address,
    squareDist=TRUE,
    8L,
    coordsV@.context_index - 1)    

range(as.matrix(D1 - D2))


#devtools::document("/home/patrick/workspace/diseasemapping/pkg/geostatsgpu")
#devtools::install("/home/patrick/workspace/diseasemapping/pkg/geostatsgpu")
#devtools::reload("/home/patrick/workspace/diseasemapping/pkg/geostatsgpu")
library(gpuR)
devtools::load_all("/home/patrick/workspace/diseasemapping/pkg/geostatsgpu")
x = as.matrix(expand.grid(1:3, 1:3))/5
coordsV = vclMatrix(x)
D3 <- vclMatrix(
    data=-0.1, 
    nrow(x), nrow(x),
    type='double'
)


maxWorkGroupSize <- gpuInfo(coordsV@.platform_index, coordsV@.device_index)$maxWorkGroupSize

.Call("cpp_gpuMatrix_custom_P", 
    coordsV@address, D3@address, 
    c(shape=1.5, range=0.2, variance = 1, nugget = 0, anisoRatio = 1, anisoAngleRadians = pi/4),
    maxWorkGroupSize, 
    coordsV@.context_index - 1)

as.matrix(D3) - as.matrix(dist(as.matrix(coordsV)))/2


.Call("stuff", coordsV@address, D3@address, as.integer(0), package='geostatspgpu')

as.matrix(D3)

range(as.matrix(D1 - D3))
