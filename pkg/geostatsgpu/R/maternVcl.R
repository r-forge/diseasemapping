#' @title Matern correlation
#'
#' @description Compute Matern covariance function
#'
#' @param x points
#' @param output a VCL matrix to store result
#' @param DofLDL a VCL vector to hold the diagonals from the LDL decomposition
#' @param param parameters
#' @param form produce variance matrix or it's cholesky
#' @return The output and DofLDL objects have been altered in-place
#' @examples
#' library(gpuR)
#' setContext(grep('gpu', listContexts()$device_type)[1])
#' x = as.matrix(expand.grid(seq(1,by=1,len=23), seq(201,by=1,len=14)))
#' coordsV = vclMatrix(x,type='float')
#' currentDevice()
#' D3 <- vclMatrix(data=-1, nrow(coordsV), nrow(coordsV),type='float')
#' 
#' myParams = c(shape=4.5, range=1.5, variance = 2, nugget = 0, anisoRatio = 2, anisoAngleRadians = pi/4)
#' 
#' system.time(v1 <- maternGpu(coordsV, D3, DofLDL = NULL, myParams,'variance', type_flag=1))
#' system.time(mmat <- geostatsp::matern(sp::SpatialPoints(as.matrix(coordsV)), myParams))
#' as.matrix(v1)[1:5,1:5]
#' mmat[1:5,1:5]
#' 
#' @export
maternGpu = function(
  x, 
  output = vclMatrix(data=0, nrow(x), nrow(x),type='float', ctx_id = x@.context_index),
  DofLDL = NULL,
  param = c(range = 1, variance = 1, shape = 1), 
  form = c("variance", "cholesky", "precision", "inverseCholesky"),
  type_flag=c(1,2)) 
{
  form = gsub("iance$|esky$|ision", "", tolower(form)[1])    
  form = c(var=1,chol=2,prec=3,inversechol=4)[form]

  if(form > 2) { warning("only form= variance or cholesky are currently implemented")}

  if(is.null(DofLDL)) {
    if(form==2) {DofLDL = gpuR::vclVector(data = -1, length=nrow(x), type='float')} 
    else {DofLDL = gpuR::vclVector(data = -1, length=1, type='float')}
  }

maxWorkGroupSize <- switch(
  deviceType(output@.device_index, output@.context_index),
  "gpu" = gpuInfo(output@.device_index, output@.context_index)$maxWorkGroupSize,
  "cpu" = cpuInfo(output@.device_index, output@.context_index)$maxWorkGroupSize,
  stop("unrecognized device type")
  )

param = geostatsp::fillParam(param) 

#file <- system.file("CL", "matern.cl", package = "geostatsgpu")
#kernel <- readChar(file, file.info(file)$size)
#junkY = junkC = gpuR::vclMatrix(matrix(1.1, 2, 2))

if (type_flag==1){
fromC = cpp_maternGpuF(
  xvt,
  output, 
  DofLDL,
  param,
  form,
  maxWorkGroupSize)

  if(form == 2) 
    {res = list(L = output, D = DofLDL, det = fromC)}
  else res = output
  res}

if (type_flag==2){
  fromC = cpp_maternGpuD(
  x,
  output, 
  DofLDL,
  param,
  form,
  maxWorkGroupSize)

  if(form == 2) 
  {res = list(L = output, D = DofLDL, det = fromC)}
  else res = output
  res}

}


