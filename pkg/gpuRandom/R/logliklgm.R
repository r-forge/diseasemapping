#' @title pre_function for likfitGpu
#'
#' @useDynLib gpuRandom
#' @export


#before start, we have SpatialPointsDataFrame and spatial model

lgmGpuObjectes <- function(modelname, mydat, type=c("double", "float")){
  
  covariates = model.matrix(modelname$model$formula, data=modelname$data)
  
  temp = model.frame(modelname$model$formula, data=modelname$data)
  
  response=temp[,as.character(attributes(terms(temp))$variables)[2]]
  
  yX=vclMatrix(cbind(response,covariates),type=type)
  
  coordsGpu<-vclMatrix(mydat@coords,type=type)
  
  n = length(response)
  p = ncol(covariates)
  
  output<-list(yX=yX, coordsGpu=coordsGpu, n=n, p=p)

  output
}

  
#' @title Estimate Log-likelihood for Gaussian random fields
#'
#' @useDynLib gpuRandom
#' @export  

  likfitGpu <- function( modelname, mydat, type=c("double", "float"),
                         paramsBatch, #vclMatrix of parameter sets,Vbatch, # matern correlation vclmatrix,diagMat, # D of cholesky decomposition
                         betas=NULL, #a vclmatrix  #given by the user or provided from formula
                         form = c("loglik", "ml", "mlFixBeta", "mlFixSigma", "reml", "remlPro"),
                         workgroupSize,
                         localSize,
                         NlocalCache,
                         verbose=FALSE){
     
     if((is.null(betas) & form == "loglik") | (is.null(betas) & form == "mlFixSigma") ) stop("need betas")
     
     yX <- lgmGpuObjectes(modelname, mydat, type)$yX
     
     coordsGpu <- lgmGpuObjectes(modelname, mydat, type)$coordsGpu   
     
     n <- lgmGpuObjectes(modelname, mydat, type)$n 
     p <- lgmGpuObjectes(modelname, mydat, type)$p  
                          
     y <-  vclMatrix(yX[,1], nrow=n, ncol=1, type = type) 
     X <-  vclMatrix(yX[,c(2:(1+p))],type = type)
                           
     form = c(loglik=1, ml=2, mlFixBeta=3, mlFixSigma=4, reml=5, remlPro=6)[form]
     
     rowbatch = nrow(paramsBatch)
     colbatch = 1     #colbatch = ncol(y)
   
     localSizechol<-localSize
     localSizechol[2]<-workgroupSize[2]
  
     
     Vbatch = vclMatrix(0, nrow(paramsBatch)*n, n, type = type)
     diagMat = vclMatrix(0, nrow(paramsBatch), n, type = type) 
     
     aTDa <- vclMatrix(0, colbatch*rowbatch, colbatch, type = gpuR::typeof(Vbatch))
     nine <- vclMatrix(0, colbatch*rowbatch, colbatch, type = gpuR::typeof(Vbatch))
     
# ########################## 1, loglik or ml(beta,sigma), given beta and sigma ##################################################

    # Vbatch=LDL^T, cholesky decomposition
    gpuRandom:::maternBatchBackend(Vbatch, coordsGpu, paramsBatch,  workgroupSize, localSize)
    gpuRandom::cholBatch(Vbatch, diagMat, numbatchD=rowbatch, Nglobal=workgroupSize, Nlocal=localSizechol, NlocalCache=NlocalCache)

    
    # temp = y-X*beta
    temp <- y - gpuRandom::gemmBatch(X, betas, 1L, 1L, colbatch, need_transpose = FALSE, workgroupSize)
    
    # L * C = temp, backsolve for C
    C <- vclMatrix(0, nrow(Vbatch), ncol(temp), type = gpuR::typeof(Vbatch))
    
    gpuRandom::backsolveBatch(C, Vbatch, temp, numbatchB=1L,  diagIsOne=TRUE,workgroupSize,  localSize,  NlocalCache)
    
    
    # one0 = C^T * D^(-1) * C = (y-X*betas)^T * V^(-1) * (y-X*betas)
    one0 <- vclMatrix(0, rowbatch*colbatch, colbatch, type = gpuR::typeof(Vbatch))
    one<- vclMatrix(0, nrow=rowbatch, ncol=colbatch, type = gpuR::typeof(Vbatch))
  
    gpuRandom:::crossprodBatchBackend(one0, C, diagMat,  invertD=TRUE,  workgroupSize, localSize, NlocalCache)

    #
    for (j in 1:colbatch){
     for(i in 1:rowbatch){
       one[i,j] <- one0[colbatch*(i-1)+j,j]
     }
    }

    #onecpu<-as.matrix(one)
    
    # part1 = n*log(sigma^2)+log |D|
   
    logD <- apply(log(diagMat),1,sum)
    
    part1 <-n*log(paramsBatch[,3]) + logD
    #replicate the part1 to do the plus operation
    part1<- matrix(part1, nrow=length(part1), ncol=colbatch, byrow=F)
    part1 <- vclMatrix(part1,type = gpuR::typeof(Vbatch))

    # n*log(sigma^2)+log |D| + one/variances
    variances<-vclMatrix(paramsBatch[,3],nrow=rowbatch, ncol=1,type = gpuR::typeof(Vbatch))
    loglik <- part1 + one/variances

   

    ###################################2, ml or ml(hatbeta,hatsigma)############################################################
    #profile = nlog hat_sigma^2 + log|D|
    # to get hat_sigma^2,  #L(a b) = (y X)
    ab <- vclMatrix(0, nrow(Vbatch), colbatch+p, type = gpuR::typeof(Vbatch))
    gpuRandom::backsolveBatch(ab, Vbatch, yX, numbatchB=1L, diagIsOne=TRUE, Nglobal=workgroupSize, Nlocal=localSize, 
                               NlocalCache)
    
    # vbatchcpu<-as.matrix(Vbatch)
    # abcpu<-as.matrix(ab)
    # yxcpu<-as.matrix(yX)
    # diagmatcpu<-as.matrix(diagMat)
    
    # temp2 = (ab)^T * D^(-1) *ab
    temp2 <- vclMatrix(0, ncol(ab)*rowbatch, ncol(ab), type = gpuR::typeof(Vbatch))
    gpuRandom:::crossprodBatchBackend(temp2, ab, diagMat, invertD=TRUE, workgroupSize, localSize, NlocalCache)
    
    #temp2cpu<-as.matrix(temp2)
    
    # b^T * D^(-1) * b = Q * P * Q^T, cholesky of a subset of temp2
    diagP <- vclMatrix(0, rowbatch, p, type = gpuR::typeof(Vbatch))
    gpuRandom:::cholBatchBackend(temp2, diagP, c(colbatch, p, colbatch, p), c(0, rowbatch, 0, p), rowbatch, workgroupSize, localSizechol, NlocalCache) 
    
    # Q * temp3 = (b^T * D^(-1) *a), backsolve for temp3
    temp3 <- vclMatrix(0, rowbatch*p, colbatch, type = gpuR::typeof(Vbatch))
    gpuRandom:::backsolveBatchBackend(temp3, temp2, temp2,
                          c(0,p,0,colbatch), c(colbatch, p, colbatch, p), c(colbatch, p, 0, colbatch),rowbatch,
                          diagIsOne=TRUE, workgroupSize, localSize, NlocalCache)
    
    # nine = temp3^T * P^(-1) * temp3,  four 2 by 2 matrices
    nine0 <- vclMatrix(0, colbatch*rowbatch, colbatch, type = gpuR::typeof(Vbatch))
    gpuRandom:::crossprodBatchBackend(nine0, temp3, diagP, invertD=TRUE,  workgroupSize, localSize, NlocalCache) ##doesn't need selecting row/col
    
   
    #extract a^TDa from temp2
    for (j in 1:colbatch){
      for (i in 1:rowbatch){
        aTDa[i,j]= temp2[(i-1)*ncol(ab)+j, j]
      }
    }
    #extract needed cells from nine0
    for (j in 1:colbatch){
      for (i in 1:rowbatch){
        nine[i,j]= nine0[(i-1)*colbatch+j, j]
      }
    }
    
    #a^TDa - nine
    two = aTDa - nine
    ml = n*log(two) + replicate(colbatch, logD)
    
    ####################################3, mlFixBeta / or ml(beta,hatsigma)##############################################
    mlFixBeta = n*log(one)+replicate(colbatch, logD)
    
    ####################################4, mlFixSigma / or ml(hatbeta,sigma)##############################################
    mlFixSigma = part1 + two/variances
    
    
    ####################################5, reml ##############################################
    #(n-p) * log sigma^2 + log |D| + log |P| + two/sigma^2
    logP <- apply(log(diagP),1,sum)
    first_part <- (n-p)*log(paramsBatch[,3]) + logD + logP
    reml <- replicate(colbatch, first_part) + two/variances
    
    
    ##################################6, remlPro or reml(hatsigma) ################################################
    #(n-p)*log two + log|D| + log|P|
    remlPro <- (n-p)*log(two) + replicate(colbatch, (logD+logP))
    
    
    if (form==1 ){
      result = loglik
    } else if (form==2){
      result =  ml
    } else if (form==3) {
      result = mlFixBeta
    } else if (form==4) {
      result= mlFixSigma
    }else if (form==5) {
      result=reml
    }else if (form==6) {
      result=remlPro
    }
    
    result
  }
    
    
    
    
    
    
    
    
  
  
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     