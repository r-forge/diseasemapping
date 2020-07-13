#' @title Estimate Log-likelihood for Gaussian random fields
#'
#' @useDynLib gpuRandom
#' @export

# likfitGpu <- function(){
#   
# }
# 
# 
# 
# 
# likfitGpu0 <- function(y, #vclmatrix of observations
#                       X, #vclmatrix of covaraites (), not #of explanatory variables?
#                       coordsGpu,
#                       paramsBatch, #vclMatrix of parameter sets
#                       Vbatch, # matern correlation vclmatrix
#                       diagMat, # D of cholesky decomposition
#                       betas, #
#                       form = c("loglik", "pro_loglik", "reml", "res_pro_loglik")
#                       minustwotimes=TRUE,
#                       workgroupSize1,
#                       workgroupSize2, #for gemmbatch
#                       localSize,
#                       NlocalCache=1000,
#                       verbose=FALSE){
# 
# 
# 
#   #1, get matern correlation matrix Vbatch=R()+mu^2*I
#   #Vbatch = vclMatrix(0, nrow(paramsBatch)*nrow(y), nrow(y), type=gpuR::typeof(paramsBatch))
#   maternBatchBackend(Vbatch, coordsGpu, paramsBatch,  Nglobal = workgroupSize1, Nlocal = localSize)
# 
#   #2, Vbatch=LDL^T, cholesky decomposition
#   #diagMat = vclMatrix(0, nrow(paramsBatch), ncol(Vbatch), type = gpuR::typeof(Vbatch))
#   cholBatchBackend(Vbatch, diagMat, Nglobal = workgroupSize1, Nlocal = localSize, NlocalCache = NlocalCache)
# 
#   # get part1 = n*log(sigma^2)+log |D|
#   n<-nrow(y)
#   part1_0 <-n*log(paramsBatch['variance']) + apply(log(diagMat),1,sum)
#   part1 <- replicate(ncol(y), part1_0)
#   
#   #3, L(a,b)=(y,X)   use backsolve function
#   # a = L^-1 * y, b = L^-1 * X
#   a <- vclMatrix(0, nrow(Vbatch), ncol(y), type = gpuR::typeof(Vbatch))  # eg, 40x2
#   b <- vclMatrix(0, nrow(Vbatch), ncol(X), type = gpuR::typeof(Vbatch))  #40x3
# 
#   backsolveBatchBackend(a, Vbatch, y, diagIsOne=TRUE, workgroupSize1, Nlocal= localSize, NlocalCache)
#   backsolveBatchBackend(b, Vbatch, X, diagIsOne=TRUE, workgroupSize1, Nlocal= localSize, NlocalCache)
# 
#   #4, D(c,e)=(a,b) backsolve
#   c <- vclMatrix(0, nrow(Vbatch), ncol(a), type = gpuR::typeof(Vbatch))   #40x2
#   e <- vclMatrix(0, nrow(Vbatch), ncol(b), type = gpuR::typeof(Vbatch))   #40x3
# 
#   backsolveBatchBackend(c, diagMat, a, diagIsOne=FALSE, workgroupSize1, Nlocal= localSize, NlocalCache)  #???
#   backsolveBatchBackend(e, diagMat, b, diagIsOne=FALSE, workgroupSize1, Nlocal= localSize, NlocalCache)  #???
# 
#   #5, get part2=(a^T*c - 2*a^T e * beta + beta^T * b^T * e * beta)/sigma^2
# 
#   aTc<- vclMatrix(0, nrow(paramsBatch), ncol(y), type = gpuR::typeof(Vbatch)) # 4 sets of 1x2 = 4x2
#   gemmBatchBackend(a,c,aTc,nrow(paramsBatch),Acolbatch=ncol(y), Bcolbatch=ncol(y), need_transpose=TRUE,Nglobal=workgroupSize2)
# 
#   aTe <- vclMatrix(0, nrow(paramsBatch), ncol(y)*ncol(e), type = gpuR::typeof(Vbatch))
#   gemmBatchBackend(a,e, aTe, nrow(paramsBatch), Acolbatch=ncol(y), Bcolbatch=1L, need_transpose=TRUE,Nglobal=workgroupSize2 )
#   
#   aTe_beta <- vclMatrix(0, nrow(paramsBatch), ncol(y), type = gpuR::typeof(Vbatch)) # 4 sets of 1x2 = 4x2
#   gemmBatchBackend(aTe,betas, aTe_beta, nrow(paramsBatch), Acolbatch=ncol(y), Bcolbatch=ncol(y), need_transpose=FALSE,Nglobal=workgroupSize2 )
#   
#   bTe <- vclMatrix(0, nrow(paramsBatch)*ncol(b), ncol(e), type = gpuR::typeof(Vbatch))  # 4 sets of 3x3
#   gemmBatchBackend(b,e, bTe, nrow(paramsBatch), Acolbatch=1L, Bcolbatch=1L, need_transpose=T,Nglobal=workgroupSize2 )
#   
#   beta_bTe<- vclMatrix(0, nrow(paramsBatch), ncol(bTe)*ncol(y), type = gpuR::typeof(Vbatch)) 
#   gemmBatchBackend(betas, bTe, beta_bTe, nrow(paramsBatch), Acolbatch=ncol(y), Bcolbatch=1L, need_transpose=T,Nglobal=workgroupSize2 )
#   
#   
#   beta_be_beta<- vclMatrix(0, nrow(paramsBatch), ncol(y), type = gpuR::typeof(Vbatch)) # 4 sets of 1x2 = 4x2
#   gemmBatchBackend(beta_bTe, betas, beta_be_beta, nrow(paramsBatch), Acolbatch=ncol(y), Bcolbatch=ncol(y), need_transpose=F, Nglobal=workgroupSize2)
#   
# 
#   #6
#   part2<- (aTc - 2*aTe_beta + beta_be_beta)/paramsBatch['variance']
# 
#   #7 
#   if (form == "loglik")
#   {result<-part1+part2}
# ##################################################################################  
#   #1 hat_betas= (b^T D^-1 b)^(-1)* b^T *c
#   
#   
#   
#   if (form == "pro_loglik"){
#     hat_betas<- 
#     
#     
#     result <- }
#   
# }
# 
# 
# 
# 





