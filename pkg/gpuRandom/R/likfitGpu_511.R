#' @title Estimate Log-likelihood for Gaussian random fields
#' @useDynLib gpuRandom
#' @export 



likfitGpu_2 <- function(modelname, 
                      data, 
                      type,
                      wholeparamsBatch, #a vclmatrix
                      betas, #a vclmatrix  #given by the user or provided from formula
                      wholeLogLik,
                      BoxCox, # an R vector
                      form = c("loglik", "ml", "mlFixSigma", "mlFixBeta", "reml", "remlPro"),
                      NparamPerIter,  # how many rows of params to be executed in each loop
                      workgroupSize,
                      localSize,
                      NlocalCache){
  
  form = c(loglik=1, ml=2, mlFixSigma=3, mlFixBeta=4, reml=5, remlPro=6)[form]
  
  covariates = model.matrix(modelname$model$formula, data=modelname$data)
  temp = model.frame(modelname$model$formula, data=modelname$data)
  y = temp[,as.character(attributes(terms(temp))$variables)[2]]
  n = length(y)
  
  coordsGpu<-vclMatrix(data@coords,type=type)
  
  
  if(BoxCox[1] != 1) {
    BoxCox = c(1, 0, setdiff(BoxCox, c(0,1)))
  }  
  
  yx = vclMatrix(cbind(y,  
                       matrix(0, n, length(BoxCox)-1),  
                       covariates), 
                 type=type)
  
  #as.matrix(yx)
  

  
  
  Ncov = ncol(covariates)
  Ndata = length(BoxCox)
  Nparam = nrow(wholeparamsBatch)
  BoxCoxGpu = vclVector(BoxCox, type=type)#  mat <- matrix(0, Ncov * Nparam,  Ndata)#  betas = vclMatrix(mat, type=type)
  variances <- vclMatrix(matrix(wholeparamsBatch[,3], nrow=Nparam, ncol=Ndata, byrow=FALSE), type=type)     # this has to be a vclMatrix not vector
  
  
  logD = vclVector(0, Nparam, type=type)
  logP = vclVector(0, Nparam, type=type)
  jacobian = vclVector(0, Ndata, type=type)
  ssqY <- vclMatrix(0, Nparam, Ndata, type=type)
  ssqbetahat <- vclMatrix(0, Nparam, Ndata, type=type)
  
  
  
  ssqYX = vclMatrix(0, ncol(yx) * NparamPerIter, ncol(yx), type=type)
  LinvYX = vclMatrix(0, nrow(yx) * NparamPerIter, ncol(yx), type=type)
  QinvSsqYx = vclMatrix(0, NparamPerIter*Ncov, Ndata, type = type)
  varMat = vclMatrix(0, nrow(yx) * NparamPerIter, nrow(coordsGpu), type=type)
  cholXVXdiag = vclMatrix(0, NparamPerIter, Ndata, type=type)
  minusTwoLogLik = vclMatrix(0, NparamPerIter, Ndata, type=type)
  
  
  gpuRandom:::likfitGpu_BackendP(
    yx,
    coordsGpu,
    paramsBatch,
    BoxCoxGpu,
    betas,
    ssqY,
    ssqbetahat,
    logD,
    logP,
    jacobian,
    NparamPerIter=10,
    workgroupSize,
    localSize,
    NlocalCache,
    verbose = 0,
    ssqYX, LinvYX,QinvSsqYx,cholXVXdiag) 
  
  
  as.vector(logD)
  as.vector(logP)
  as.vector(jacobian)
  as.matrix(ssqbetahat) 
  as.matrix(ssqY)[1:20,]
  as.matrix(ssqYX)
  
  # resid^T V^(-1) resid, resid = Y - X betahat = two
  two <- ssqY - ssqbetahat
  
  if(form ==2){ #ml  
    # = n*log(two/n) + logD + jacobian +n + n*log(2*pi) 
    
    temp1 = n*log(two/n)
    matrix_vector_sumBackend(temp1, logD, jacobian, n + n*log(2*pi), minusTwoLogLik, workgroupSize)
    wholeLogLik = -0.5*minusTwoLogLik
    
  }else if(form==4){ #mlFixBeta
    # two/variances + n*log(variances) + logD + jacobian + n*log(2*pi)
    temp1 = two/variances + n*log(variances)
    
    matrix_vector_sumBackend(temp1, logD, jacobian, n*log(2*pi), minusTwoLogLik, workgroupSize)
    
    wholeLogLik = -0.5 * minusTwoLogLik
    
  }else if(form==5){ #reml
    #first_part <- (n-p)*log(variances) + logD + logP  
    #minusTwoLogLik=first_part + two/variances +jacobian + n*log(2*pi)
    temp1 = two/variances + (n-p)*log(variances)
    temp2 = logD + logP  
    gpuRandom:::matrix_vector_sumBackend(temp1, temp2, jacobian, n*log(2*pi), minusTwoLogLik, workgroupSize)
    wholeLogLik = -0.5 * minusTwoLogLik
    
  }else if(form==6){ # remlPro
    # minusTwoLogLik= (n-p)*log(two/(n-p)) + logD + logP + jacobian + n*log(2*pi) + n-p
    temp1 = (n-p)*log(two/(n-p))
    temp2 = logD + logP 
    gpuRandom:::matrix_vector_sumBackend(temp1, temp2, jacobian, n-p+ n*log(2*pi), minusTwoLogLik, workgroupSize)
    
    wholeLogLik = -0.5 * minusTwoLogLik
  } 
  
  
  
  
}




