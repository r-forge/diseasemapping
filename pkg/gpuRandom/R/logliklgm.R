#' @title pre_function1 for likfitGpu
#'
#' @useDynLib gpuRandom
#' @export


#before start, we have spatial model and SpatialPointsDataFrame 


lgmGpuObjectes1 <- function(modelname, mydat, type=c("double", "float")){
    
    covariates = model.matrix(modelname$model$formula, data=modelname$data)
    temp = model.frame(modelname$model$formula, data=modelname$data)
    response=temp[,as.character(attributes(terms(temp))$variables)[2]]
    n = length(response)
    p = ncol(covariates)
    yX=vclMatrix(cbind(response,covariates),type=type)   # y <-vclMatrix(yX[,1], nrow=n, ncol=1, type = type) # X <-vclMatrix(yX[,c(2:(1+p))],type = type)
    coordsGpu<-vclMatrix(mydat@coords,type=type)
    
    
    output<-list(yX=yX, coordsGpu=coordsGpu, n=n, p=p)
    
    output
}

#' @title pre_function2 for likfitGpu
#'
#' @useDynLib gpuRandom
#' @export

lgmGpuObjectes2 <- function(rowbatch, n, p, type=c("double", "float")){
    
    colbatch = 1
    
    Vbatch = vclMatrix(0, rowbatch*n, n, type = type)
    diagMat = vclMatrix(0, rowbatch, n, type = type) 
    logD<- vclVector(0, length=rowbatch,type=type)
    ab <- vclMatrix(0, nrow(Vbatch), 1+p, type = type)
    A <- vclMatrix(0, nrow(Vbatch), colbatch, type = type)
    temp0 <- vclMatrix(0, rowbatch, colbatch, type = type)   
    temp1 <- vclMatrix(0, rowbatch, colbatch, type = type)   
    temp2 <- vclMatrix(0, (1+p)*rowbatch, (1+p), type = type)
    diagP <- vclMatrix(0, rowbatch, p, type = type)
    temp3 <- vclMatrix(0, rowbatch*p, colbatch, type = type)
    nine0 <- vclMatrix(0, colbatch*rowbatch, colbatch, type = type)
    aTDa <- vclMatrix(0, colbatch*rowbatch, colbatch, type = type)
    ssqBeta <- vclMatrix(0, nrow=rowbatch, ncol=1, type=type)
    logP <- vclVector(0, length=rowbatch, type=type)
    
    
    output<-list(Vbatch=Vbatch, diagMat=diagMat, logD=logD, ab=ab, A=A, temp0=temp0, temp1=temp1, temp2=temp2, diagP=diagP, 
                 temp3=temp3, nine0=nine0, aTDa=aTDa, ssqBeta=ssqBeta, logP=logP)
    
    output
}




#' @title Estimate Log-likelihood for Gaussian random fields middle step
#'
#' @useDynLib gpuRandom
#' @export  

likfitGpu_0 <- function(yX, n, p, coordsGpu, 
                        type=c("double", "float"),
                        Vbatch,
                        diagMat,
                        logD,
                        ab,
                        A,
                        temp0,
                        temp1,
                        temp2,
                        diagP,
                        temp3,
                        nine0,
                        aTDa,
                        ssqBeta,
                        logP,
                        paramsBatch, 
                        betas=NULL, #a vclmatrix  #given by the user or provided from formula
                        jacobian,
                        form = c("loglik", "ml", "mlFixSigma", "mlFixBeta", "reml", "remlPro"),# minustwotimes=TRUE,
                        workgroupSize,
                        localSize,
                        NlocalCache){
    
    if((is.null(betas) & form == "loglik") | (is.null(betas) & form == "mlFixSigma") ) stop("need betas")
    
    form = c(loglik=1, ml=2, mlFixSigma=3, mlFixBeta=4, reml=5, remlPro=6)[form]
    
    rowbatch = nrow(paramsBatch)
    colbatch = 1     #colbatch = ncol(y)
    localSizechol<-localSize
    localSizechol[2]<-workgroupSize[2]
    
    ########################### 1, loglik or ml(beta,sigma), given beta and sigma #############################
    #get matern matrix
    gpuRandom::maternBatch(Vbatch, coordsGpu, paramsBatch,  workgroupSize, localSize)
    #Vbatch=LDL^T, cholesky decomposition
    gpuRandom::cholBatch(Vbatch, diagMat, numbatchD=rowbatch, Nglobal=workgroupSize, Nlocal=localSizechol, NlocalCache=NlocalCache)
    #logD <- apply(log(diagMat),1,sum)
    gpuRandom:::rowsumBackend(diagMat, logD, type="row", log=1)    
    variances<-vclMatrix(paramsBatch[,3],nrow=rowbatch, ncol=1,type = type)   
    #L(a b) = (y X)
    gpuRandom::backsolveBatch(ab, Vbatch, yX, numbatchB=1L, diagIsOne=TRUE, Nglobal=workgroupSize, Nlocal=localSize, NlocalCache)
    # temp2 = (ab)^T * D^(-1) *ab
    gpuRandom::crossprodBatch(temp2, ab, diagMat, invertD=TRUE, workgroupSize, localSize, NlocalCache)    
    
    #extract a^TDa from temp2
    for (j in 1:colbatch){
        for (i in 1:rowbatch){
            aTDa[i,j]= temp2[(i-1)*ncol(ab)+j, j]
        }
    }
    
    if(form == 1 | form == 3){       
        # a^TD^(-1)b * beta = temp0
        gemmBatch2backend(temp2, betas, temp0,
                          c(0,0,0), # transposeA, transposeB, transposeC
                          c(0,1,(1+p),1,p,(1+p)), c(0,p,p,0,1,1), #select from temp2 and betas
                          c(0,1,1,0,1,1), #select from output
                          c(rowbatch,1,0,0,1,0), # batches: nRow, nCol, recycleArow, recycleAcol, recycleB row col
                          c(workgroupSize,localSize), #workgroupSize  global 0 1 2, local 0 1 2  
                          c(NlocalCache,NlocalCache), # cacheSizeA, cacheSizeB,                
                          FALSE )
        
        
        # b * beta = temp1 (n x 1), to get ssqBeta
        gemmBatch2backend(ab, betas, temp1,
                          c(0,0,0), # transposeA, transposeB, transposeC
                          c(0,n,n,1,p,(1+p)), c(0,p,p,0,1,1), #select from ab and betas
                          c(0,n,n,0,1,1), #select from output
                          c(rowbatch,1,0,0,1,0), # batches: nRow, nCol, recycleArow, recycleAcol, recycleB row col
                          c(workgroupSize,localSize), #workgroupSize  global 0 1 2, local 0 1 2  
                          c(NlocalCache,NlocalCache), # cacheSizeA, cacheSizeB,                
                          FALSE )
        
        # ssqBeta = temp1^T D^(-1) temp1 = beta T*bT D^(-1) b*beta
        gpuRandom::crossprodBatch(ssqBeta, temp1, diagMat,  invertD=TRUE,  workgroupSize, localSize, NlocalCache)
        
        # to get one, see the notes
        one0 = aTDa - 2*temp0 + ssqBeta
    }
    # to get hat_sigma^2 when beta is not given by users
    # b^T * D^(-1) * b = Q * P * Q^T, cholesky of a subset (right bottom) of temp2
    gpuRandom:::cholBatchBackend(temp2, diagP, c(colbatch, p, colbatch, p), c(0, rowbatch, 0, p), rowbatch, workgroupSize, localSizechol, NlocalCache)     
    # Q * temp3 = (b^T * D^(-1) *a), backsolve for temp3    2 by 1
    gpuRandom:::backsolveBatchBackend(temp3, temp2, temp2, c(0,p,0,colbatch), c(colbatch, p, colbatch, p), c(colbatch, p, 0, colbatch),
                                      rowbatch,diagIsOne=TRUE, workgroupSize, localSize, NlocalCache)     
    # nine0 = temp3^T * P^(-1) * temp3,  four 1 by 1 matrices
    gpuRandom::crossprodBatch(nine0, temp3, diagP, invertD=TRUE,  workgroupSize, localSize, NlocalCache) ##doesn't need selecting row/col
    
    two = aTDa - nine0  
    
    if(form == 1 | form == 4){     
        part1 <- n*log(variances) + logD    #part1 = n*log(sigma^2)+log |D|
    }  
    gpuRandom:::rowsumBackend(diagP, logP, type="row",log=1)  #logP <- apply(log(diagP),1,sum)
    
    
    
    
    if(form == 1 ){ #loglik
        # n*log(sigma^2) + log |D| + one/variances # result <- part1 + one0/variances + n*log(2*pi)
        Result = list(minusTwoLogLik=(part1 + one0/variances + n*log(2*pi))*jacobian, 
                      LogLik = -0.5*(part1 + one0/variances + n*log(2*pi))*jacobian,
                      ssqBeta=ssqBeta, ssqX=NULL, ssqY=aTDa, logD=logD, logP=logP)      
        
    }else if(form == 2) {#ml  result = n*log(two) +logD + n*log(2*pi) + n
        Result = list(minusTwoLogLik=(n*log(two/n) +logD + n*log(2*pi) + n)*jacobian,
                      LogLik = -0.5*(n*log(two/n) +logD + n*log(2*pi) + n)*jacobian,
                      ssqBeta=0, ssqX=nine0, ssqY=aTDa, logD=logD, logP=logP)           
        
    }else if(form == 3){ # mlFixSigma/ or ml(beta,hatsigma)result = n*log(one0/n)+logD + n*log(2*pi) + n
        Result = list(minusTwoLogLik=(n*log(one0/n)+logD + n*log(2*pi) + n)*jacobian, 
                      LogLik = -0.5*(n*log(one0/n)+logD + n*log(2*pi) + n)*jacobian,
                      ssqBeta=ssqBeta, ssqX=nine0, ssqY=aTDa, logD=logD, logP=logP)   
        
    }else if(form == 4){ # mlFixBeta / or ml(hatbeta,sigma) result = part1 + two/variances + n*log(2*pi) 
        Result = list(minusTwoLogLik=(part1 + two/variances + n*log(2*pi))*jacobian, 
                      LogLik = -0.5*(part1 + two/variances + n*log(2*pi))*jacobian,
                      ssqBeta=0, ssqX=nine0, ssqY=aTDa, logD=logD, logP=logP)            
        
    }else if(form == 5){ #reml
        first_part <- (n-p)*log(variances) + logD + logP  #result <- first_part + two/variances + n*log(2*pi) 
        Result = list(minusTwoLogLik=(first_part + two/variances + n*log(2*pi))*jacobian, 
                      LogLik = -0.5*(first_part + two/variances + n*log(2*pi))*jacobian,
                      ssqBeta=0, ssqX=nine0, ssqY=aTDa, logD=logD, logP=logP)            
        
    }else if(form == 6){ #remlPro  #(n-p)*log two + log|D| + log|P|, result <- (n-p)*log(two/(n-p)) + logD+logP + n*log(2*pi) + n-p
        Result = list(minusTwoLogLik=((n-p)*log(two/(n-p)) + logD+logP + n*log(2*pi) + n-p)*jacobian, 
                      LogLik = -0.5*((n-p)*log(two/(n-p)) + logD+logP + n*log(2*pi) + n-p)*jacobian,
                      ssqBeta=0, ssqX=nine0, ssqY=aTDa, logD=logD, logP=logP)                  
    }
    
    Result
    
}


#' @title Estimate Log-likelihood for Gaussian random fields
#'
#' @useDynLib gpuRandom
#' @export 
likfitGpu <- function(modelname, mydat, type=c("double", "float"), 
                      bigparamsBatch, 
                      betas=NULL, #a vclmatrix  #given by the user or provided from formula
                      BoxCox,
                      form = c("loglik", "ml", "mlFixSigma", "mlFixBeta", "reml", "remlPro"),# minustwotimes=TRUE,
                      groupsize,  # how many rows of params to be executed each loop
                      workgroupSize,
                      localSize,
                      NlocalCache,
                      verbose=FALSE){
    
    output1 <- lgmGpuObjectes1(modelname, mydat, type=type)
    
    output2 <- lgmGpuObjectes2(groupsize, output1$n, output1$p, type=type)
    
    totalnumbersets <- nrow(bigparamsBatch)
    index <- c(1:groupsize)
    
    LogLik <- vclVector(rep(0, totalnumbersets),type=type)
    # box cox transform
    
    yX <- output1$yX
    
    jacobian = 0
    
    if(BoxCox != 1) {
        
        if(abs(BoxCox - 1 ) < 0.001) {
            jacobian=0  # BoxCox close to 1, don't transform
        }
        else{# box cox is not one.
            jacobian = -2*(BoxCox-1)* sum(log(yX[,1]))  
            
            if(is.nan(jacobian))
                warning("boxcox shouldnt be used with negative data")
            
            if(abs(BoxCox)<0.001) {
                yX[,1] = log(yX[,1]) 
            }else if(abs(BoxCox-1)>0.001) {
                yX[,1] <- ((yX[,1]^BoxCox) - 1)/BoxCox
            }
        } 
    } # end have box cox
    
    
    
    for(i in 0:(totalnumbersets/groupsize) ){    # can do >8990 sets of parameters       
        
        if(groupsize*(i+1) < totalnumbersets +1){
            
            paramsBatchcpu<-bigparamsBatch[index + i*groupsize,]
            
            paramsBatch <- vclMatrix(paramsBatchcpu, type=type)
            
            resulti <- likfitGpu_0(yX,
                                   output1$n, 
                                   output1$p, 
                                   output1$coordsGpu, 
                                   type=type,
                                   output2$Vbatch,
                                   output2$diagMat,
                                   output2$logD,
                                   output2$ab,
                                   output2$A,
                                   output2$temp0,
                                   output2$temp1,
                                   output2$temp2,
                                   output2$diagP,
                                   output2$temp3,
                                   output2$nine0,
                                   output2$aTDa,
                                   output2$ssqBeta,
                                   output2$logP,
                                   paramsBatch, 
                                   betas= betas, #a vclmatrix  #given by the user or provided from formula
                                   jacobian,
                                   form = form,# minustwotimes=TRUE,
                                   workgroupSize,
                                   localSize,
                                   NlocalCache)
            
            
            replace(LogLik,index + i*groupsize, resulti$LogLik)
            
            
            
        }   
    }
    
    LogLik
}




