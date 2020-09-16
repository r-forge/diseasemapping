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
  likfitGpu0 <- function(y, #vclmatrix of observations
                         X, #vclmatrix of covaraites (), not #of explanatory variables?
                         coordsGpu,
                         paramsBatch, #vclMatrix of parameter sets
                         Vbatch, # matern correlation vclmatrix
                         diagMat, # D of cholesky decomposition
                         betas, # vclmatricesform = c("loglik", "pro_loglik", "reml", "pro_reml")
                         NlocalCache,
                         workgroupSize,
                         localSize,
                         verbose=FALSE){



     rowbatch = nrow(paramsBatch)
     colbatch = ncol(y)
     n = nrow(y)
     p = ncol(X)

# ##########################loglik########################################################

    #1, Vbatch=LDL^T, cholesky decomposition
    gpuRandom:::maternBatchBackend(Vbatch, coordsGpu, paramsBatch,  workgroupSize, localSize)
    Vbatch<-gpuRandom::cholBatch(Vbatch, diagMat, numbatchD=rowbatch, Nglobal=workgroupSize, Nlocal=localSize, NlocalCache=NlocalCache)$L
    diagMat<-gpuRandom::cholBatch(Vbatch, diagMat, numbatchD=rowbatch, Nglobal=workgroupSize, Nlocal=localSize, NlocalCache=NlocalCache)$diag
    

    #2, temp = y-X*beta
    temp <- y - gpuRandom::gemmBatch(X, betas, rowbatch, 1L, colbatch, need_transpose = FALSE, workgroupSize)

    #3, C = L^(-1) * temp,   L * C = temp, backsolve
    C <- vclMatrix(0, nrow(Vbatch), ncol(temp), type = gpuR::typeof(Vbatch))

    gpuRandom:::backsolveBatchBackend(C, Vbatch, temp,
                          c(0,n,0,ncol(y)),   c(0,n,0,n),  c(0,n,0,ncol(y)),
                          rowbatch,  diagIsOne=TRUE,
                          workgroupSize,  localSize,  NlocalCache)

    #4, result0 = C^T * D^(-1) * C = (y-X*betas)^T * V^(-1) * (y-X*betas)
    result0 <- vclMatrix(0, rowbatch*colbatch, colbatch, type = gpuR::typeof(Vbatch))

    gpuRandom:::crossprodBatchBackend(result0, C, diagMat,  invertD=T,  workgroupSize, localSize, NlocalCache)

    part2<- vclMatrix(0, nrow=rowbatch, ncol=colbatch, type = gpuR::typeof(Vbatch))

    #need edit
    for (j in 1:colbatch){
     for(i in 1:rowbatch){
      part2[i,j] <- result0[colbatch*(i-1)+j,j]
     }
    }

    
    
    #5, part1 = n*log(sigma^2)+log |D|
   
    logD <- apply(log(diagMat),1,sum)
    
    part1_0 <-n*log(paramsBatch[,3]) + logD
    #replicate the part1 to do the plus operation
    part1_0<- matrix(part1_0, nrow=length(part1_0), ncol=colbatch, byrow=F)
    part1 <- vclMatrix(part1_0,type = gpuR::typeof(Vbatch))

    #6,
    variances<-vclMatrix(paramsBatch[,3],nrow=rowbatch, ncol=1,type = gpuR::typeof(Vbatch))
    loglik <- part1 + part2/variances


    loglik

  }
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
   #1, get matern correlation matrix Vbatch=R()+mu^2*I
    #Vbatch = vclMatrix(0, nrow(paramsBatch)*nrow(y), nrow(y), type=gpuR::typeof(paramsBatch))
   # maternBatchBackend(Vbatch, coordsGpu, paramsBatch,  Nglobal = workgroupSize1, Nlocal = localSize)

   #2, Vbatch=LDL^T, cholesky decomposition
    #diagMat = vclMatrix(0, nrow(paramsBatch), ncol(Vbatch), type = gpuR::typeof(Vbatch))
   # cholBatchBackend(Vbatch, diagMat, Nglobal = workgroupSize1, Nlocal = localSize, NlocalCache = NlocalCache)
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
#   if (form == "loglik")
#   { part2<- (aTc - 2*aTe_beta + beta_be_beta)/paramsBatch['variance']
#
#     result<- part1+part2}
# ##################################################################################
#   #1 hat_betas= (b^T D^-1 b)^(-1)* b^T *c = b^(-1) D c
#   bTc<- gpuRandom::gemmBatch(b, c, nrow(paramsBatch), 1L, ncol(y), TRUE, workgroupSize2)
#
#   # Df=b  backsolve f=D^(-1)*b
#
#
#   bDb<-vclMatrix(0, nrow(paramsBatch)*ncol(b), ncol(b), type = gpuR::typeof(Vbatch))
#   crossprodBatchBackend(bDb, b, diagMat, invertD=T, workgroupSize1, Nlocal, NlocalCache)
#
#
#
#
#   hat_betas <-gemmBatch(bDbinverse, bTc, nrow(paramsBatch), 1L, 1L, F, workgroupSize2)
#
#   #2, n*log hat_sigma^2 + log |D|
#   # hat_sigma^2= cTDTc-aTc
#   if (form == "pro_loglik"){
#
#
#
#     result <- }
#
 









