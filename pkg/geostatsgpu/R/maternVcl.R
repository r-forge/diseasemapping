#' @title Matern correlation
#'
#' @description Compute Matern covariance function
#'
#' @param x points
#' @param param parameters
#' @param output VCL matrix to store result
#' @param type produce variance matrix or it's cholesky
#' @param DofLDL a VCL vector to hold the diagonals from the LDL decomposition
#' @return The output and DofLDL objects have been altered in-place
#' @examples
#' x = as.matrix(expand.grid(seq(1,by=1,len=23), seq(201,by=1,len=14)))
#' coordsV = gpuR::vclMatrix(VCLmatrix(x))
#' D3 <- vclMatrix(
#'     data=-1, 
#'     nrow(coordsV), nrow(coordsV),
#'     type=typeof(coordsV)
#' )
#' 
#' myParams = c(shape=4.5, range=1.5, variance = 2, nugget = 0, anisoRatio = 2, anisoAngleRadians = pi/4)
#' 
#' system.time(v1 <- maternGpu(coordsV, myParams, D3, 'variance'))
#' system.time(mmat <- geostatsp::matern(sp::SpatialPoints(as.matrix(coordsV)), myParams))
#' as.matrix(v1)[1:5,1:5]
#' mmat[1:5,1:5]
#' 
#' @export
maternGpu = function(
  x, 
  param = c(range = 1, variance = 1, shape = 1), 
  output,
  type = c("variance", "cholesky", "precision", "inverseCholesky"),
  DofLDL = NULL
  ) {

  type = gsub("iance$|esky$|ision", "", tolower(type)[1])    
  type = c(var=1,chol=2,prec=3,inversechol=4)[type]

  if(type > 2) { 
    warning("only type= variance or cholesky are currently implemented")
  }


if(is.data.frame(x)) {
  x = as.matrix(x)
}
if(is.matrix(x)) {
  x = vclMatrix(x)
}

if(missing(output)) {
  output = vclMatrix(
    data=0, 
    nrow(x), nrow(x),
    type='double',
    ctx_id = x@.context_index)
}

 if(is.null(DofLDL)) {
    if(type==2) {
      DofLDL = gpuR::vclVector(data = -1, length=nrow(x), type='double', ctx_id = x@.context_index)
    } else {
      DofLDL = gpuR::vclVector(data = -1, length=1, type='double', ctx_id = x@.context_index)
    }
 }

maxWorkGroupSize <- switch(
  deviceType(output@.platform_index, output@.device_index),
  "gpu" = gpuInfo(output@.platform_index, output@.device_index)$maxWorkGroupSize,
  "cpu" = cpuInfo(output@.platform_index, output@.device_index)$maxWorkGroupSize,
  stop("unrecognized device type")
  )

param = geostatsp::fillParam(param)

#file <- system.file("CL", "matern.cl", package = "geostatsgpu")
#kernel <- readChar(file, file.info(file)$size)


#junkY = junkC = gpuR::vclMatrix(matrix(1.1, 2, 2))

fromC = geostatsgpu:::cpp_maternGpu(
  output, 
  DofLDL,
  x,
  param[c('range','shape','variance','nugget','anisoRatio','anisoAngleRadians')],
  type,
  maxWorkGroupSize)

  if(type == 2) {
    res = list(L = output, D = DofLDL, logDet = fromC)
  }
  else res = output

  res
}


