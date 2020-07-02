#' @title Estimate Log-likelihood for Gaussian random fields
#'
#' @useDynLib gpuRandom
#' @export









# # 1, log |sigma^2 * V|=n*log (sigma^2) + log |D|
# loglikGpu = function(data,  # An object of class SpatialPointsDataFrame
#                      trend,
#                      paramsBatch, # a vclMatrix of parameters shape, range, variance, nugget, anisoRatio, anisoAngleRadians
#                      minustwotimes=TRUE,
#                      workgroupSize,
#                      localSize,
#                      NlocalCache=1000, 
#                      verbose=FALSE
#                      ){
#   
#   #coordinates
#   
#   coordsGpu = vclMatrix(data@coords, nrow(data@coords), ncol(data@coords), type=gpuR::typeof(paramBatch))
#   
#   # get matern correlation matrix V=R+mu^2*I
#   outputBatchV = vclMatrix(0, nrow(paramBatch)*nrow(coordsGpu), nrow(coordsGpu), type=gpuR::typeof(paramsBatch))
#   
#   maternBatchBackend(outputBatchV, coordsGpu, paramBatch,  Nglobal = workgroupSize, Nlocal = localSize)
#   
#   # do cholesky factorization of V
#   
#   outputChol = vclMatrix(0, nrow(outputBatchV), ncol(outputBatchV), type = gpuR::typeof(outputBatchV))
#   
#   diagMat = vclMatrix(0, nrow(paramBatch), ncol(outputBatchV), type = gpuR::typeof(outputBatchV))
#   
#  
#   cholBatchBackend(outputChol, diagMat, Nglobal = workgroupSize,  Nlocal = localSize, NlocalCache = NlocalCache)
#   
#   
#   # part 1 of loglikelihood
#   n<-nrow(outputBatchV)
#   diagM<-as.matrix(diagMat)
#   part1<-n*log(paramsBatch['variance']) + apply(log(diagM),1,sum)   # log(gpumatrix)
#   
#   
#   # 2, (y-X\beta)^\T V^(-1)(y-X\beta) / sigma^2= a^Tc-2 beta b^T c + beta^T b^T e beta /sigma^2
#  
#    # 2(1) L(a b)=(y X)
#   observations = all.vars(trend)[1]
#   covariates = model.frame(trend, data.frame(data))
#   
#   y = covariates[,	observations]
#   y<-vclMatrix(y)
#   covariates = model.matrix(trend, covariates, drop.unused.levels = FALSE)
#   X<-vclMatrix(covariates)
#   
#   
#   solvea=vclMatrix(0, nrow(paramBatch)*nrow(outputChol),1,type=gpuR::typeof(outputChol))
#   backsolveBatchBackend(solvea, outputChol, y,  diagIsOne = TRUE,   Nglobal = NglobalChol,Nlocal = NlocalChol, NlocalCache = NlocalCache)
#   
#   solveb=vclMatrix(0,nrow(paramBatch)*nrow(outputChol), ncol(X),type=gpuR::typeof(outputChol))
#   backsolveBatchBackend(solveb, outputChol, y,  diagIsOne = TRUE,   Nglobal = NglobalChol,Nlocal = NlocalChol, NlocalCache = NlocalCache)
#   
#   
#   #2(2)  c=D^-1 a   e=D^-1 b
#   
#   c<-vclMatrix(0, nrow(paramBatch)*nrow(outputChol),1,type=gpuR::typeof(outputChol))
#   e<-vclMatrix(0, nrow(paramBatch)*nrow(outputChol), ncol(X),type=gpuR::typeof(outputChol))
#   
#   multiplyDiagonalBatchBackend(c, diagMat, solvea,inverse=TRUE,Nglobal,Nlocal)
#   multiplyDiagonalBatchBackend(e, diagMat, solveb,inverse=TRUE,Nglobal, Nlocal) 
#   
#   #2(3) a^T D^-1 a - 2beta^T b^T c + beta^T b^T e beta
#   
#   first<-vclMatrix(0, nrow(paramBatch),1,type=gpuR::typeof(outputChol))
#   crossprodBatchBackend(first,solvea,diagMat,invertD=TRUE, Nglobal,Nlocal,  NlocalCache)
#   
#   second<-
#   
# }















