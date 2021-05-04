#' @title pre_function for likfitGpu
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







likfitGpu_0 <- function(Vbatch,
                        coordsGpu, 
                        paramsBatchperIter, 
                        diagMat, 
                        logD_temp,   #viennacl::vector
                        logD_i,
                        ab,
                        yX, #y1,y2,y3,X            
                        temp2,
                        aTDa,   # ssqY
                        betas, #a vclmatrix  //given by the user or provided from formula, default=null
                        temp00,
                        temp0,
                        temp1,
                        ssqBeta0,
                        ssqBeta, # 
                        diagP,
                        temp3,
                        ssqbetahat0,
                        ssqbetahat,      # ssqbetahat
                        logP_temp,  # viennacl::vector
                        logP_i,
                        Qinverse,
                        identity,
                        QPQinverse,
                        betahat,  # p*rowbatch   colbatch                            
                        variances,   # must be rowbatch * colbatch matrix !!!
                        form_temp,
                        form_temp1,
                        jacobian,    #a vclmatrix 
                        LogLik,
                        n, 
                        p, 
                        rowbatch, 
                        colbatch,   
                        form, # c(loglik=1, ml=2, mlFixSigma=3, mlFixBeta=4, reml=5, remlPro=6)
                        workgroupSize, localSize, localSizechol, NlocalCache){
  
  
  #get matern matrix
  gpuRandom::maternBatch(Vbatch, coordsGpu, paramsBatchperIter,  workgroupSize, localSize)
  #Vbatch=LDL^T, cholesky decomposition
  gpuRandom::cholBatch(Vbatch, diagMat, numbatchD=rowbatch, Nglobal=workgroupSize, Nlocal=localSizechol, NlocalCache=NlocalCache)
  
  #logD_temp <- apply(log(diagMat),1,sum)   half log determinant of V
  gpuRandom:::rowsumBackend(diagMat, logD_temp, type="row", log=TRUE)   
  logD_i[] <- logD_temp
  
  #L(a1,a2,a3, b) = (y1,y2,y3, X)
  gpuRandom::backsolveBatch(ab, Vbatch, yX, numbatchB=1L, diagIsOne=TRUE, Nglobal=workgroupSize, Nlocal=localSize, NlocalCache)
  
  # temp2 = (ab)^T * D^(-1) *ab  = (Y X)^T V^{-1} (YX)
  gpuRandom::crossprodBatch(temp2, ab, diagMat, invertD=TRUE, workgroupSize, localSize, NlocalCache)    
  
  
  
  # SSQY
  #// SSQY = a^TDa from temp2       1x3
  for (j in 1:colbatch){
    for (i in 1:rowbatch){
      aTDa[i,j]= temp2[(i-1)*ncol(ab)+j, j]
    }
  }
  
  
  # beta is supplied 
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
  # X^T V^(-1) X = b^T * D^(-1) * b = Q * P * Q^T, cholesky of a subset (right bottom) of temp2
  gpuRandom:::cholBatchBackend(temp2, diagP, c(colbatch, p, colbatch, p), c(0, rowbatch, 0, p), rowbatch, workgroupSize, localSizechol, NlocalCache)     
  
  
  # store determinant for REML  
  gpuRandom:::rowsumBackend(diagP, logP_temp, type="row",log=1L)  
  logP_i[] <- logP_temp
  
  
  # ssqbetahat = (X betahat)^T V^(-1) X betahat 
  #  Q * temp3 = (b^T * D^(-1) *a), backsolve for temp3    2 by 1
  gpuRandom:::backsolveBatchBackend(temp3, temp2, temp2, 
                                    c(0,p,0,colbatch), 
                                    c(colbatch, p, colbatch, p), 
                                    c(colbatch, p, 0, colbatch),
                                    rowbatch,diagIsOne=TRUE, workgroupSize, localSize, NlocalCache)   
  
  
  # ssqbetahat0 = temp3^T * P^(-1) * temp3,  four 1 by 1 matrices      3x3
  gpuRandom::crossprodBatch(ssqbetahat0, temp3, diagP, invertD=TRUE,  workgroupSize, localSize, NlocalCache) 
  
  
  # TO DO (long list), save ssqbetahat matrix, row batches and column batches
  # save diagonals of ssqbetahat, for standard errors of beta hat
  # extract needed cells from ssqbetahat0       1x3
  for (j in 1:colbatch){
    for (i in 1:rowbatch){
      ssqbetahat[i,j]= ssqbetahat0[(i-1)*colbatch+j, j]
    }
  }   
  
  #// resid^T V^(-1) resid, resid = Y - X betahat 
  two = aTDa - ssqbetahat   #// log ssqResid   
  
  
  # betahat = (X^T V^(-1) X)^(-1)  X^T V^(-1) Y
  # (X^T V^(-1) X) =Q P Q^T 
  # calculate betahat = Q^(-T) P^(-1) Q^(-1) * b^T D^(-1)a
  # first: get Q^(-1)    pxp
  gpuRandom:::backsolveBatchBackend(Qinverse, temp2, identity, 
                                    c(0, p, 0, p), 
                                    c(colbatch, p, colbatch, p), 
                                    c(0, p, 0, p),
                                    numbatchB=1L, diagIsOne=TRUE, workgroupSize, localSize, NlocalCache)   
  
  
  
  #// Q^(-T) P^(-1) Q^(-1) = QPQinverse 
  gpuRandom::crossprodBatch(QPQinverse, Qinverse, diagP,  invertD=TRUE,  workgroupSize, localSize, NlocalCache) 
  # must be a batch of square matrices 
  
  
  #// QPQinverse * b^T D^(-1)a = betahat
  gemmBatch2backend(QPQinverse, temp2, betahat,
                    c(0,0,0), # transposeA, transposeB, transposeC
                    c(0,p,p,0,p,p), 
                    c(colbatch, p,(colbatch+p), 0, colbatch, (colbatch+p)), #select from temp2 and betas
                    c(0,p,p,0,colbatch,colbatch), #select from output
                    c(rowbatch,1,0,0,0,0), # batches: nRow, nCol, recycleArow, recycleAcol, recycleBrow, col
                    c(workgroupSize,localSize), #workgroupSize  global 0 1 2, local 0 1 2  
                    c(NlocalCache,NlocalCache), # cacheSizeA, cacheSizeB,                
                    FALSE )
  
  
  
  if(form == 1 | form == 4){     # variances are provided
    # get form_temp <- n*log(variances) + logD_temp    , form_temp = n*log(sigma^2)+log |D|
    new_temp = n*log(variances)
    gpuRandom:::matrix_vector_sumBackend(new_temp, logD_temp, form_temp, byrow = TRUE, workgroupSize)
    
  } 
  
  if (form ==5 | form ==6){
    logD_plusP = logD_temp + logP_temp
  }
  
  
  
  
  
  if (form==1){ #loglik
    # form_temp + one/variances +jacobian + n*log(2*pi)
    minusTwoLogLik = form_temp + one/variances + jacobian + n*log(2*pi)
    
    LogLik[,] = -0.5*minusTwoLogLik
    
  }else if(form ==2){ #ml  
    # = n*log(two) +logD + n*log(2*pi) + n
    new_temp = n*log(two/n)
    gpuRandom:::matrix_vector_sumBackend(new_temp, logD_temp, form_temp1, byrow = TRUE, workgroupSize)
    minusTwoLogLik= form_temp1 + jacobian + n*log(2*pi) + n 
    LogLik[,] = -0.5*minusTwoLogLik
    
  }else if(form==3){ # mlFixSigma
    # n*log(one/n)+ logD +jacobian + n*log(2*pi) + n
    new_temp = n*log(one/n)
    gpuRandom:::matrix_vector_sumBackend(new_temp, logD_temp, form_temp1, byrow = TRUE, workgroupSize)  
    minusTwoLogLik= form_temp1 + jacobian + n*log(2*pi) + n 
    LogLik[,] = -0.5*minusTwoLogLik
    
  }else if(form==4){ #mlFixBeta
    #form_temp + two/variances +jacobian + n*log(2*pi)
    
    minusTwoLogLik= form_temp + two/variances +jacobian + n*log(2*pi)
    LogLik[,] = -0.5 * minusTwoLogLik
    
  }else if(form==5){ #reml
    #first_part <- (n-p)*log(variances) + logD + logP  
    #minusTwoLogLik=first_part + two/variances +jacobian + n*log(2*pi)
    
    new_temp = (n-p)*log(variances)
    gpuRandom:::matrix_vector_sumBackend(new_temp, logD_plusP, form_temp1, byrow = TRUE, workgroupSize)
    
    minusTwoLogLik= form_temp1 + two/variances +jacobian + n*log(2*pi)
    LogLik[,] = -0.5 * minusTwoLogLik
    
  }else if(form==6){ # remlPro
    # minusTwoLogLik= (n-p)*log(two/(n-p)) + logD + logP + jacobian + n*log(2*pi) + n-p
    new_temp = (n-p)*log(two/(n-p))
    gpuRandom:::matrix_vector_sumBackend(new_temp, logD_plusP, form_temp1, byrow = TRUE, workgroupSize)
    
    minusTwoLogLik= form_temp1 +jacobian + n*log(2*pi) + n-p 
    LogLik[,] = -0.5 * minusTwoLogLik
  } 
  
  
}











#' @title Estimate Log-likelihood for Gaussian random fields
#' @useDynLib gpuRandom
#' @export 
likfitGpu <- function(modelname, mydat, type=c("double", "float"), 
                      completeparamsBatch, #a vclmatrix
                      betas=NULL, #a vclmatrix  #given by the user or provided from formula
                      logD,
                      logP,
                      BoxCox, # an R vector
                      form = c("loglik", "ml", "mlFixSigma", "mlFixBeta", "reml", "remlPro"),
                      groupsize,  # how many rows of params to be executed in each loop
                      workgroupSize,
                      localSize,
                      NlocalCache){
  
  form = c(loglik=1, ml=2, mlFixSigma=3, mlFixBeta=4, reml=5, remlPro=6)[form]
  localSizechol<-localSize
  localSizechol[2]<-workgroupSize[2]
  totalnumbersets <- nrow(completeparamsBatch)     
  ncols <- ncol(completeparamsBatch)
  
  output1 <- lgmGpuObjectes1(modelname, mydat, type=type)
  colbatch<- as.integer(length(BoxCox)+1) 
  n<-output1$n
  p<-output1$p
  
  ssqBeta <- vclMatrix(0, groupsize, colbatch, type=type)
  ssqbetahat  <- vclMatrix(0, groupsize, colbatch, type = type)
  aTDa <- vclMatrix(0, groupsize, colbatch, type = type)
  logD_temp <- vclVector(0, length=groupsize, type=type)
  logP_temp <- vclVector(0, length=groupsize, type=type)
  betahat <- vclMatrix(0, groupsize*p, colbatch, type = type)
  
  
  Vbatch <- vclMatrix(0, groupsize*n, n, type = type)
  diagMat <- vclMatrix(0, groupsize, n, type = type) 
  ab <- vclMatrix(0, groupsize*n, colbatch+p, type = type)
  temp2 <- vclMatrix(0, (colbatch+p)*groupsize, (colbatch+p), type = type)
  temp00 <- vclMatrix(0, groupsize*colbatch, colbatch, type = type)   
  temp0 <- vclMatrix(0, groupsize, colbatch, type = type)   
  temp1 <- vclMatrix(0, groupsize*n, colbatch, type = type)   
  ssqBeta0 <- vclMatrix(0, groupsize*colbatch, colbatch, type=type)
  one <- vclMatrix(0, groupsize, colbatch, type=type)
  diagP <- vclMatrix(0, groupsize, p, type = type)
  temp3 <- vclMatrix(0, groupsize*p, colbatch, type = type)
  ssqbetahat0 <- vclMatrix(0, colbatch*groupsize, colbatch, type = type)
  two  <- vclMatrix(0, groupsize, colbatch, type = type)
  Qinverse  <- vclMatrix(0, groupsize*p, p, type = type)
  identity <- gpuR::identity_matrix(p, type = type)
  QPQinverse  <- vclMatrix(0, groupsize*p, p, type = type)
  form_temp <- vclMatrix(0, groupsize, colbatch, type = type)
  form_temp1 <- vclMatrix(0, groupsize, colbatch, type = type)
  finalLogLik <- vclMatrix(0, totalnumbersets, colbatch, type=type)
  
  yXcpu <- output1$yXcpu
  
  # get jacobian matrix
  jacobian = -2*(BoxCox-1)* sum(log(yXcpu[,1])) 
  closetooneindex <- which(abs(BoxCox - 1 ) < 0.001)
  jacobian[closetooneindex] = 0    
  jacobian <- c(jacobian, 0)   
  jacobian<- vclMatrix(matrix(jacobian, nrow=groupsize, ncol=length(jacobian), byrow=TRUE), type=type) # make it from a vector to a matrix!!!
  
  
  # box cox transform   
  transformed_y = matrix(0,n,length(BoxCox))    
  for (i in 1:length(BoxCox)){
    transformed_y[ ,i] <- ((yXcpu[ ,1]^BoxCox[i]) - 1)/BoxCox[i]
  }      
  closetozeroindex <- which(abs(BoxCox)<0.001)
  transformed_y[ ,closetozeroindex] = log(yXcpu[,1])  
  transformed_y[ ,closetooneindex] = yXcpu[,1]   
  yX <- vclMatrix(cbind(transformed_y,yXcpu),type=type)
  
  
  variances <- vclMatrix(matrix(completeparamsBatch[,3], nrow=totalnumbersets, ncol=colbatch, byrow=FALSE), type=type)     # this has to be a vclMatrix not vector
  
  
  
  N <- totalnumbersets/groupsize
  
  for(i in 1:N ){    # can do >8990 sets of parameters       
    
    paramsBatch_i <-gpuR::block(completeparamsBatch, rowStart = as.integer(1L + (i-1)*groupsize), rowEnd = as.integer(i * groupsize), colStart = 1L, colEnd = ncols)
    params_foruse <- deepcopy(paramsBatch_i)
    
    logD_i <- gpuR::slice(logD, start=as.integer(1L + (i-1)*groupsize), end=as.integer(i * groupsize))
    logP_i <- gpuR::slice(logP, start=as.integer(1L + (i-1)*groupsize), end=as.integer(i * groupsize))
    
    LogLik_i <-gpuR::block(finalLogLik, rowStart = as.integer(1L + (i-1)*groupsize), rowEnd = as.integer(i * groupsize), colStart = 1L, colEnd = colbatch)
    
    variances_i <- gpuR::block(variances, rowStart = as.integer(1L + (i-1)*groupsize), rowEnd = as.integer(i * groupsize), colStart = 1L, colEnd = colbatch)
    
    likfitGpu_0(Vbatch,
                output1$coordsGpu,
                params_foruse, 
                diagMat, 
                logD_temp,   #viennacl::vector
                logD_i,
                ab,
                yX, #y1,y2,y3,X            
                temp2,
                aTDa,   # ssqY
                betas=betas, #a vclmatrix  //given by the user or provided from formula, default=null
                temp00,
                temp0,
                temp1,
                ssqBeta0,
                ssqBeta, # 
                diagP,
                temp3,
                ssqbetahat0,
                ssqbetahat,      # ssqbetahat
                logP_temp,  # viennacl::vector
                logP_i,
                Qinverse,
                identity,
                QPQinverse,
                betahat,  # p*rowbatch   colbatch                            
                variances_i,   # must be rowbatch * colbatch matrix !!!
                form_temp,
                form_temp1,
                jacobian,    #a vclmatrix 
                LogLik_i,
                n, 
                p, 
                groupsize, 
                colbatch,   
                form=form, # c(loglik=1, ml=2, mlFixSigma=3, mlFixBeta=4, reml=5, remlPro=6)
                workgroupSize, localSize, localSizechol, NlocalCache)
    
    
  }
  
  finalLogLik
}















