#' @title Estimate Log-likelihood for Gaussian random fields
#'
#' @useDynLib gpuRandom
#' @export


loglikGpu <- function(y, #vclmatrix of observations
                      X, #vclmatrix of covaraites (), not #of explanatory variables?
                      coordsGpu,
                      paramsBatch, #vclMatrix of parameter sets
                      betas, #???
                      minustwotimes=TRUE,
                      workgroupSize1,
                      workgroupSize2,
                      localSize,
                      NlocalCache=1000, 
                      verbose=FALSE){
  
  
  
  #1, get matern correlation matrix Vbatch=R()+mu^2*I
  Vbatch = vclMatrix(0, nrow(paramsBatch)*nrow(y), nrow(y), type=gpuR::typeof(paramsBatch))
  
  maternBatchBackend(Vbatch, coordsGpu, paramsBatch,  Nglobal = workgroupSize1, Nlocal = localSize)
  
  #2, Vbatch=LDL^T, cholesky decomposition
  diagMat = vclMatrix(0, nrow(paramsBatch), ncol(Vbatch), type = gpuR::typeof(Vbatch))
  
  cholBatchBackend(Vbatch, diagMat, Nglobal = workgroupSize1, Nlocal = localSize, NlocalCache = NlocalCache)
  
  # get part1 = n*log(sigma^2)+log |D|
  n<-nrow(y)
  part1<-n*log(paramsBatch['variance']) + apply(log(diagMat),1,sum) 
  
  #3, L(a,b)=(y,X)   use backsolve function
  # a = L^-1 y, b = L^-1 X
  a <- vclMatrix(0, nrow(Vbatch), ncol(y), type = gpuR::typeof(Vbatch))  # eg, 40x2
  b <- vclMatrix(0, nrow(Vbatch), ncol(X), type = gpuR::typeof(Vbatch))  #40x3
  
  backsolveBatchBackend(a, Vbatch, y, diagIsOne=TRUE, workgroupSize1, Nlocal= localSize, NlocalCache)
  backsolveBatchBackend(b, Vbatch, X, diagIsOne=TRUE, workgroupSize1, Nlocal= localSize, NlocalCache)
  
  # to do, ab = backsoolveBactchBackend(ab, Vbatch, cbind(y,X), ...)

  #4, D(c,e)=(a,b) backsolve
  c <- vclMatrix(0, nrow(Vbatch), ncol(a), type = gpuR::typeof(Vbatch))   #40x2
  e <- vclMatrix(0, nrow(Vbatch), ncol(b), type = gpuR::typeof(Vbatch))   #40x3

  # to do ce
  
  backsolveBatchBackend(c, diagMat, a, diagIsOne=FALSE, workgroupSize1, Nlocal= localSize, NlocalCache)  #???
  backsolveBatchBackend(e, diagMat, b, diagIsOne=FALSE, workgroupSize1, Nlocal= localSize, NlocalCache)  #???
  
  #5, get part2=(a^T*c - 2*a^T e * beta + beta^T * b^T * e * beta)/sigma^2

  #aTce = [ aT c, aT e]
  gemmBatchBackend(
    aTce
    ab, 1, 1, # ab[,1:ncol(y)]
    ce, 1, ncol(ce), 
    1 # transpose
    )


  aTe = gemmBatchBackend(
    ab, 1, ncol(y),# ab[,1:ncol(y)]
    ce, ncol(y)+1, ncol(y) + ncol(X),
    1)

  ac<- vclMatrix(0, ncol(y)*nrow(paramsBatch), ncol(y), type = gpuR::typeof(Vbatch)) # 4 sets of 2x2 =8x2
  ae<- vclMatrix(0, ncol(y)*nrow(paramsBatch), ncol(b), type = gpuR::typeof(Vbatch)) # 4 sets of 2x3 =8x3
  be<- vclMatrix(0, ncol(X)*nrow(paramsBatch), ncol(b), type = gpuR::typeof(Vbatch)) # 4 sets of 3x3 =12x3
  
  gemmBatchBackend(a,c,ac,nrow(paramsBatch),need_transpose=1,Nglobal=workgroupSize2)
  gemmBatchBackend(a,e,ae,nrow(paramsBatch),need_transpose=1,Nglobal=workgroupSize2)
  gemmBatchBackend(b,e,be,nrow(paramsBatch),need_transpose=1,Nglobal=workgroupSize2)
  
     #beta 4 sets of 3x2
  ae_beta<- vclMatrix(0, ncol(y)*nrow(paramsBatch), ncol(beta), type = gpuR::typeof(Vbatch))  # 4 sets of 2x2
  gemmBatchBackend(ae, beta, ae_beta, nrow(paramsBatch), need_transpose=0, Nglobal=workgroupSize2)
  
  beta_be<- vclMatrix(0, nrow(paramsBatch)*ncol(beta), ncol(be), type = gpuR::typeof(Vbatch)) # 4 sets of 2x3
  gemmBatchBackend(beta, be, beta_be, nrow(paramsBatch), need_transpose=1, Nglobal=workgroupSize2)
  
  
  beta_be_beta<- vclMatrix(0, nrow(paramsBatch)*ncol(beta), ncol(beta), type = gpuR::typeof(Vbatch)) # 4 sets of 2x2
  gemmBatchBackend(beta_be, beta, beta_be_beta, nrow(paramsBatch), need_transpose=0, Nglobal=workgroupSize2)
  
  #6,
  part2<- ac - 2*ae_beta
  
  
}







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















