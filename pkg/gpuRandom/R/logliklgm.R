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
    yXcpu=cbind(response,covariates)   # y <-vclMatrix(yX[,1], nrow=n, ncol=1, type = type) # X <-vclMatrix(yX[,c(2:(1+p))],type = type)
    coordsGpu<-vclMatrix(mydat@coords,type=type)
    
    
    output<-list(yXcpu=yXcpu, coordsGpu=coordsGpu, n=n, p=p)
    
    output
}

#' @title pre_function2 for likfitGpu
#'
#' @useDynLib gpuRandom
#' @export

lgmGpuObjectes2 <- function(rowbatch, colbatch, n, p, type=c("double", "float")){
    
    Vbatch = vclMatrix(0, rowbatch*n, n, type = type)
    diagMat = vclMatrix(0, rowbatch, n, type = type) 
    #logD<- vclVector(0, length=rowbatch,type=type)
    ab <- vclMatrix(0, rowbatch*n, colbatch+p, type = type)
    temp00 <- vclMatrix(0, rowbatch*colbatch, colbatch, type = type)   
    temp0 <- vclMatrix(0, rowbatch, colbatch, type = type)   
    temp1 <- vclMatrix(0, rowbatch*n, colbatch, type = type)   
    temp2 <- vclMatrix(0, (colbatch+p)*rowbatch, (colbatch+p), type = type)
    diagP <- vclMatrix(0, rowbatch, p, type = type)
    temp3 <- vclMatrix(0, rowbatch*p, colbatch, type = type)
    nine0 <- vclMatrix(0, colbatch*rowbatch, colbatch, type = type)
    nine <- vclMatrix(0, rowbatch, colbatch, type = type)
    aTDa <- vclMatrix(0, rowbatch, colbatch, type = type)
    ssqBeta0 <- vclMatrix(0, rowbatch*colbatch, colbatch, type=type)
    ssqBeta <- vclMatrix(0, rowbatch, colbatch, type=type)
    #logP <- vclVector(0, length=rowbatch, type=type)
    
    
    output<-list(Vbatch=Vbatch, diagMat=diagMat, #logD=logD, 
                 ab=ab, temp00=temp00, temp0=temp0,
                 temp1=temp1, temp2=temp2, diagP=diagP, temp3=temp3, 
                 nine0=nine0, nine=nine, aTDa=aTDa, ssqBeta0=ssqBeta0, ssqBeta=ssqBeta) #logP=logP)
    
    output
}




#' @title Estimate Log-likelihood for Gaussian random fields middle step
#'
#' @useDynLib gpuRandom
#' @export  

likfitGpu_0 <- function(yX,  #y1,y2,y3,X
                        n, p, coordsGpu, 
                        rowbatch, colbatch,
                        type=c("double", "float"),
                        Vbatch,
                        diagMat, #logD,
                        ab,
                        temp00,
                        temp0,
                        temp1,
                        temp2,
                        diagP,
                        temp3,
                        nine0,
                        nine,
                        aTDa,
                        ssqBeta0,
                        ssqBeta, #logP,
                        paramsBatch, 
                        betas=NULL, #a vclmatrix  #given by the user or provided from formula
                        jacobian,
                        form = c("loglik", "ml", "mlFixSigma", "mlFixBeta", "reml", "remlPro"),# minustwotimes=TRUE,
                        workgroupSize,
                        localSize,
                        NlocalCache){
    
    if((is.null(betas) & form == "loglik") | (is.null(betas) & form == "mlFixSigma") ) stop("need betas")
    
    form = c(loglik=1, ml=2, mlFixSigma=3, mlFixBeta=4, reml=5, remlPro=6)[form]
    
    # rowbatch = nrow(paramsBatch)
    # colbatch = ncol(yX)-1
    localSizechol<-localSize
    localSizechol[2]<-workgroupSize[2]
    
    ########################### 1, loglik or ml(beta,sigma), given beta and sigma #############################
    #get matern matrix
    gpuRandom::maternBatch(Vbatch, coordsGpu, paramsBatch,  workgroupSize, localSize)
    #Vbatch=LDL^T, cholesky decomposition
    gpuRandom::cholBatch(Vbatch, diagMat, numbatchD=rowbatch, Nglobal=workgroupSize, Nlocal=localSizechol, NlocalCache=NlocalCache)
    #logD <- apply(log(diagMat),1,sum)
    # gpuRandom:::rowsumBackend(diagMat, logD, type="row", log=1)   
    logD <- apply(log(diagMat),1,sum)
    logD <- vclMatrix(matrix(logD, nrow=rowbatch, ncol=colbatch, byrow=FALSE), type=type)
    variances <- vclMatrix(matrix(paramsBatch[,3], nrow=rowbatch, ncol=colbatch, byrow=FALSE), type=type)     # this has to be a vclMatrix not vector
    
    #L(a1,a2,a3, b) = (y1,y2,y3, X)
    gpuRandom::backsolveBatch(ab, Vbatch, yX, numbatchB=1L, diagIsOne=TRUE, Nglobal=workgroupSize, Nlocal=localSize, NlocalCache)
    # temp2 = (ab)^T * D^(-1) *ab
    gpuRandom::crossprodBatch(temp2, ab, diagMat, invertD=TRUE, workgroupSize, localSize, NlocalCache)    
    
    #extract a^TDa from temp2       1x3
    for (j in 1:colbatch){
        for (i in 1:rowbatch){
            aTDa[i,j]= temp2[(i-1)*ncol(ab)+j, j]
        }
    }
    
    if(form == 1 | form == 3){       
        # a^TD^(-1)b * beta = temp0
        gemmBatch2backend(temp2, betas, temp00,
                          c(0,0,0), # transposeA, transposeB, transposeC
                          c(0,colbatch,(colbatch+p),colbatch,p,(colbatch+p)), c(0,p,p,0,colbatch,colbatch), #select from temp2 and betas
                          c(0,colbatch,colbatch,0,colbatch,colbatch), #select from output
                          c(rowbatch,1,0,0,1,0), # batches: nRow, nCol, recycleArow, recycleAcol, recycleBrow, col
                          c(workgroupSize,localSize), #workgroupSize  global 0 1 2, local 0 1 2  
                          c(NlocalCache,NlocalCache), # cacheSizeA, cacheSizeB,                
                          FALSE )
        # extract a^TD^(-1)b * beta from temp00      1x3
        for (j in 1:colbatch){
            for (i in 1:rowbatch){
                temp0[i,j]= temp00[(i-1)*colbatch+j, j]
            }
        }
        
        
        # b * beta = temp1 (n x colbatch), to get ssqBeta
        gemmBatch2backend(ab, betas, temp1,
                          c(0,0,0), # transposeA, transposeB, transposeC
                          c(0,n,n,colbatch,p,(colbatch+p)), c(0,p,p,0,colbatch,colbatch), #select from ab and betas
                          c(0,n,n,0,colbatch,colbatch), #select from output
                          c(rowbatch,1,0,0,1,0), # batches: nRow, nCol, recycleArow, recycleAcol, recycleB row col
                          c(workgroupSize,localSize), #workgroupSize  global 0 1 2, local 0 1 2  
                          c(NlocalCache,NlocalCache), # cacheSizeA, cacheSizeB,                
                          FALSE )
        
        # ssqBeta = temp1^T D^(-1) temp1 = beta T*bT D^(-1) b*beta    3x3
        gpuRandom::crossprodBatch(ssqBeta0, temp1, diagMat,  invertD=TRUE,  workgroupSize, localSize, NlocalCache)
        
        # extract beta T*bT D^(-1) b*beta from temp00       1x3
        for (j in 1:colbatch){
            for (i in 1:rowbatch){
                ssqBeta[i,j]= ssqBeta0[(i-1)*colbatch+j, j]
            }
        }  
        
        
        # to get one, see the notes 
        one = aTDa - 2*temp0 + ssqBeta      # 1x3
    }
    
    # to get hat_sigma^2 when beta is not given by users
    # b^T * D^(-1) * b = Q * P * Q^T, cholesky of a subset (right bottom) of temp2
    gpuRandom:::cholBatchBackend(temp2, diagP, 
                                 c(colbatch, p, colbatch, p), 
                                 c(0, rowbatch, 0, p), 
                                 rowbatch, workgroupSize, localSizechol, NlocalCache)     
    # Q * temp3 = (b^T * D^(-1) *a), backsolve for temp3    2 by 1
    gpuRandom:::backsolveBatchBackend(temp3, temp2, temp2, 
                                      c(0,p,0,colbatch), 
                                      c(colbatch, p, colbatch, p), 
                                      c(colbatch, p, 0, colbatch),
                                      rowbatch,diagIsOne=TRUE, workgroupSize, localSize, NlocalCache)     
    # nine0 = temp3^T * P^(-1) * temp3,  four 1 by 1 matrices      3x3
    gpuRandom::crossprodBatch(nine0, temp3, diagP, invertD=TRUE,  workgroupSize, localSize, NlocalCache) 
    
    # extract needed cells from nine0       1x3
    for (j in 1:colbatch){
        for (i in 1:rowbatch){
            nine[i,j]= nine0[(i-1)*colbatch+j, j]
        }
    }   
    
    two = aTDa - nine 
    
    
    
    if(form == 1 | form == 4){     
        part1 <- n*log(variances) + logD    #part1 = n*log(sigma^2)+log |D|
    }  
    
    #gpuRandom:::rowsumBackend(diagP, logP, type="row",log=1)  #
    logP <- apply(log(diagP),1,sum)
    logP <- vclMatrix(matrix(logP, nrow=length(logP), ncol=colbatch, byrow=FALSE), type=type)
    
    
    if(form == 1 ){ #loglik
        # n*log(sigma^2) + log |D| + one/variances # result <- part1 + one/variances + n*log(2*pi)
        Result = list(minusTwoLogLik= part1 + one/variances + n*log(2*pi)+jacobian, 
                      LogLik = -0.5*(part1 + one/variances + n*log(2*pi)+jacobian),
                      ssqBeta=ssqBeta, ssqX=NULL, ssqY=aTDa, logD=logD, logP=logP)      
        
    }else if(form == 2) {#ml  result = n*log(two) +logD + n*log(2*pi) + n
        Result = list(minusTwoLogLik= n*log(two/n) +logD + n*log(2*pi) + n +jacobian,
                      LogLik = -0.5*(n*log(two/n) +logD + n*log(2*pi) + n + jacobian),
                      ssqBeta=0, ssqX=nine, ssqY=aTDa, logD=logD, logP=logP)           
        
    }else if(form == 3){ # mlFixSigma/ or ml(beta,hatsigma)result = n*log(one/n)+logD + n*log(2*pi) + n
        Result = list(minusTwoLogLik=n*log(one/n)+logD + n*log(2*pi) + n+jacobian, 
                      LogLik = -0.5*(n*log(one/n)+logD + n*log(2*pi) + n+jacobian),
                      ssqBeta=ssqBeta, ssqX=nine, ssqY=aTDa, logD=logD, logP=logP)   
        
    }else if(form == 4){ # mlFixBeta / or ml(hatbeta,sigma) result = part1 + two/variances + n*log(2*pi) 
        Result = list(minusTwoLogLik= part1 + two/variances + n*log(2*pi)+jacobian, 
                      LogLik = -0.5*(part1 + two/variances + n*log(2*pi)+jacobian),
                      ssqBeta=0, ssqX=nine, ssqY=aTDa, logD=logD, logP=logP)            
        
    }else if(form == 5){ #reml
        first_part <- (n-p)*log(variances) + logD + logP  #result <- first_part + two/variances + n*log(2*pi) 
        Result = list(minusTwoLogLik=first_part + two/variances + n*log(2*pi)+jacobian, 
                      LogLik = -0.5*(first_part + two/variances + n*log(2*pi)+jacobian),
                      ssqBeta=0, ssqX=nine, ssqY=aTDa, logD=logD, logP=logP)            
        
    }else if(form == 6){ #remlPro  #(n-p)*log two + log|D| + log|P|, result <- (n-p)*log(two/(n-p)) + logD+logP + n*log(2*pi) + n-p
        Result = list(minusTwoLogLik= (n-p)*log(two/(n-p)) + logD + logP + n*log(2*pi) + n-p + jacobian, 
                      LogLik = -0.5*((n-p)*log(two/(n-p)) + logD + logP + n*log(2*pi) + n-p + jacobian),
                      ssqBeta=0, ssqX=nine, ssqY=aTDa, logD=logD, logP=logP)                  
    }
    
    Result
    
}


#' @title Estimate Log-likelihood for Gaussian random fields
#' @useDynLib gpuRandom
#' @export 
likfitGpu <- function(modelname, mydat, type=c("double", "float"), 
                      bigparamsBatchcpu, 
                      betas=NULL, #a vclmatrix  #given by the user or provided from formula
                      BoxCox, # an R vector
                      form = c("loglik", "ml", "mlFixSigma", "mlFixBeta", "reml", "remlPro"),# minustwotimes=TRUE,
                      groupsize,  # how many rows of params to be executed each loop
                      workgroupSize,
                      localSize,
                      NlocalCache,
                      verbose=FALSE){
    
    output1 <- lgmGpuObjectes1(modelname, mydat, type=type)
    colbatch<- length(BoxCox)+1   
    #lgmGpuObjectes2 <- function(rowbatch, colbatch, n, p, type=c("double", "float")){  
    output2 <- lgmGpuObjectes2(groupsize, colbatch, output1$n, output1$p, type=type)
    
    yXcpu <- output1$yXcpu
    
    # box cox transform   
    jacobian = -2*(BoxCox-1)* sum(log(yXcpu[,1])) 
    closetooneindex <- which(abs(BoxCox - 1 ) < 0.001)
    jacobian[closetooneindex] = 0    
    jacobian <- c(jacobian, 0)   
    if(is.nan(jacobian))
        warning("boxcox shouldnt be used with negative data")
    jacobian<- matrix(jacobian, nrow=groupsize, ncol=length(jacobian), byrow=TRUE)
    
    
    transformed_y = matrix(0,output1$n,length(BoxCox))    
    for (i in 1:length(BoxCox)){
        transformed_y[ ,i] <- ((yXcpu[ ,1]^BoxCox[i]) - 1)/BoxCox[i]
    }      
    closetozeroindex <- which(abs(BoxCox)<0.001)
    transformed_y[ ,closetozeroindex] = log(yXcpu[,1])     
    yX <- vclMatrix(cbind(transformed_y,yXcpu),type=type)
    
    
    totalnumbersets <- nrow(bigparamsBatchcpu)      
    #LogLik <- vclVector(rep(0, totalnumbersets),type=type)
    #LogLik <- vclMatrix(0, totalnumbersets, colbatch, type=type)
    loopindex <- c(1:groupsize)  
    
    paramsBatchcpu<-bigparamsBatchcpu[loopindex,]
    paramsBatch <- vclMatrix(paramsBatchcpu, type=type)
    
    
    result0 <- likfitGpu_0(yX,
                           output1$n, 
                           output1$p, 
                           output1$coordsGpu, 
                           groupsize, colbatch,
                           type=type,
                           output2$Vbatch,
                           output2$diagMat, #output2$logD,
                           output2$ab,
                           output2$temp00,
                           output2$temp0,
                           output2$temp1,
                           output2$temp2,
                           output2$diagP,
                           output2$temp3,
                           output2$nine0,
                           output2$nine,
                           output2$aTDa,
                           output2$ssqBeta0,
                           output2$ssqBeta, #output2$logP,
                           paramsBatch, 
                           betas= betas, #a vclmatrix  #given by the user or provided from formula
                           jacobian,
                           form = form,# minustwotimes=TRUE,
                           workgroupSize,
                           localSize,
                           NlocalCache)
    
    
    LogLik_Result <- result0$LogLik
    
    
    
    
    for(i in 1:(totalnumbersets/groupsize) ){    # can do >8990 sets of parameters       
        
        if(groupsize*(i+1) < totalnumbersets +1){
            
            paramsBatchcpu<-bigparamsBatchcpu[loopindex + i*groupsize,]
            
            paramsBatch <- vclMatrix(paramsBatchcpu, type=type)
            
            
            resulti <- likfitGpu_0(yX,
                                   output1$n, 
                                   output1$p, 
                                   output1$coordsGpu, 
                                   groupsize, colbatch,
                                   type=type,
                                   output2$Vbatch,
                                   output2$diagMat, #output2$logD,
                                   output2$ab,
                                   output2$temp00,
                                   output2$temp0,
                                   output2$temp1,
                                   output2$temp2,
                                   output2$diagP,
                                   output2$temp3,
                                   output2$nine0,
                                   output2$nine,
                                   output2$aTDa,
                                   output2$ssqBeta0,
                                   output2$ssqBeta, #output2$logP,
                                   paramsBatch, 
                                   betas= betas, #a vclmatrix  #given by the user or provided from formula
                                   jacobian,
                                   form = form,# minustwotimes=TRUE,
                                   workgroupSize,
                                   localSize,
                                   NlocalCache)
            
            
            LogLik_Result <- rbind(LogLik_Result, resulti$LogLik)
            
            
        }   
    }
    
    LogLik_Result
}



