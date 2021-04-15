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




#' @title Estimate Log-likelihood for Gaussian random fields
#' @useDynLib gpuRandom
#' @export 
likfitGpu <- function(modelname, mydat, type=c("double", "float"), 
                      bigparamsBatch, #a vclmatrix 
                      betas=NULL, #a vclmatrix  #given by the user or provided from formula
                      BoxCox, # an R vector
                      form = c("loglik", "ml", "mlFixSigma", "mlFixBeta", "reml", "remlPro"),# minustwotimes=TRUE,
                      groupsize,  # how many rows of params to be executed each loop
                      workgroupSize,
                      localSize,
                      NlocalCache,
                      verbose=FALSE){
  
  
  form = c(loglik=1, ml=2, mlFixSigma=3, mlFixBeta=4, reml=5, remlPro=6)[form]
  
  output1 <- lgmGpuObjectes1(modelname, mydat, type=type)
  colbatch<- length(BoxCox)+1   
  
  yXcpu <- output1$yXcpu
  
  # get jacobian matrix
  jacobian = -2*(BoxCox-1)* sum(log(yXcpu[,1])) 
  closetooneindex <- which(abs(BoxCox - 1 ) < 0.001)
  jacobian[closetooneindex] = 0    
  jacobian <- c(jacobian, 0)   
  jacobian<- vclMatrix(matrix(jacobian, nrow=groupsize, ncol=length(jacobian), byrow=TRUE), type=type) # make it from a vector to a matrix!!!
  
  
  # box cox transform   
  transformed_y = matrix(0,output1$n,length(BoxCox))    
  for (i in 1:length(BoxCox)){
    transformed_y[ ,i] <- ((yXcpu[ ,1]^BoxCox[i]) - 1)/BoxCox[i]
  }      
  closetozeroindex <- which(abs(BoxCox)<0.001)
  transformed_y[ ,closetozeroindex] = log(yXcpu[,1])     
  yX <- vclMatrix(cbind(transformed_y,yXcpu),type=type)
  
  totalparams <- nrow(bigparamsBatch)
  
  bigvariances <- vclMatrix(matrix(bigparamsBatch[,3], nrow=totalparams, ncol=colbatch, byrow=FALSE), type=type)  
  ssqBetaR <- vclMatrix(0, totalparams, colbatch, type=type)
  ssqXR <- vclMatrix(0, totalparams, colbatch, type=type)
  ssqYR <- vclMatrix(0, totalparams, colbatch, type=type)
  logDR <- vclVector(0, totalparams, type=type)
  logPR <- vclVector(0, totalparams, type=type)
  finalbetahatR <- vclMatrix(0, totalparams*output1$p, colbatch, type=type)
  finalLogLik <- vclMatrix(0, totalparams, colbatch, type=type)
  
  
  likfitGpu_Backend(output1$coordsGpu,
                    bigparamsBatch,
                    yX,       # a vclMatrix of n * colbatch
                    betas,
                    bigvariances,  # a vclMatrix of nrow(bigparamsBatch) * colbatch
                    jacobian,      # a vclMatrix of groupsize * colbatch
                    ssqBetaR,      # a vclMatrix
                    ssqXR,        # a vclMatrix
                    ssqYR,        # a vclMatrix
                    logDR,         # a vclVector
                    logPR,         # a vclVector
                    finalbetahatR,      # a vclMatrix
                    finalLogLik,   # a vclMatrix of nrow(bigparamsBatch) * colbatch
                    output1$n, 
                    output1$p, 
                    groupsize,
                    colbatch,
                    form,
                    workgroupSize, 
                    localSize, 
                    NlocalCache)
  
  list(ssqX=ssqXR, ssqY=ssqYR, logD = logDR, logP = logPR, LogLik = finalLogLik)
  
  
}



