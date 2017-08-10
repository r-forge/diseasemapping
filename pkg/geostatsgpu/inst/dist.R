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
  devtools::load_all("/home/patrick/workspace/diseasemapping/pkg/geostatsgpu")
  library(gpuR)
  x = as.matrix(expand.grid(1:2, 1:2))/5
  coordsV = gpuR::vclMatrix(x)
  D3 <- vclMatrix(
      data=-0.1, 
      
      nrow(x), nrow(x),
      type='double'
  )
  
  distMat = as.matrix(dist(as.matrix(coordsV)))
  
  maxWorkGroupSize <- gpuInfo(coordsV@.platform_index, coordsV@.device_index)$maxWorkGroupSize
  
  myParams = c(shape=1.5, range=1, variance = 1, nugget = 0, anisoRatio = 1, anisoAngleRadians = pi/4)
  
  (nuround = round(myParams['shape']+0.5))
  (mu = myParams['shape'] - nuround)
  unlist(.C('Rtemme_gamma', mu, 0.1, 0.1, 0.1, 0.1))
  
  .Call("cpp_gpuMatrix_custom_P", 
      coordsV@address, D3@address, 
      myParams,
      maxWorkGroupSize, 
      coordsV@.context_index - 1)

  as.matrix(D3)#[1:3, ]
  
  matrix(.C("bessel_K_p", 
      as.integer(nuround),  as.double(mu),  as.double(myParams['shape']),
      as.double(distMat), as.double(distMat), 
      as.integer(length(distMat)), PACKAGE='geostatsgpu')[[5]], ncol(distMat)) 
  
  
besselK(distMat[1:3,1:9], myParams['shape'])

as.matrix(geostatsp::matern(sp::SpatialPoints(as.matrix(coordsV)), myParams))[1:3, 1:9]


as.matrix(D3) / geostatsp::matern(
    sp::SpatialPoints(as.matrix(coordsV)), myParams
    )


.Call("stuff", coordsV@address, D3@address, as.integer(0), package='geostatspgpu')

as.matrix(D3)

range(as.matrix(D1 - D3))
