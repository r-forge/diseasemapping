#' @title Estimate Log-likelihood for Gaussian random fields
#' @useDynLib gpuRandom
#' @export


#before start, we have a spatial model 
likfitGpu_2 <- function(spatialmodel,     #data,
                        type = c("double","float"),
                        ParamsGpu, #a vclmatrix, consists of all the parameters
                        betas=NULL, #a vclmatrix  #given by the user or provided from formula
                        BoxCox, # an R vector, will always be c(1,0,.....)
                        form = c("loglik", "ml", "mlFixSigma", "mlFixBeta", "reml", "remlPro"),
                        NparamPerIter,  # how many sets of params to be evaluated in each loop
                        minustwotimes=TRUE,
                        Nglobal,
                        Nlocal,
                        NlocalCache){
  
  form = c(loglik=1, ml=2, mlFixSigma=3, mlFixBeta=4, reml=5, remlPro=6)[form]
  
  covariates = model.matrix(spatialmodel$model$formula, data=spatialmodel$data)
  temp = model.frame(spatialmodel$model$formula, data=spatialmodel$data)
  y = temp[,as.character(attributes(terms(temp))$variables)[2]]
  n = length(y)
  variances <- vclVector(ParamsGpu[,3],type=type)     # a vclvector
  
  
  
  if(BoxCox[1] != 1) {
    BoxCox = c(1, 0, setdiff(BoxCox, c(0,1)))
  }
  
  Ncov = ncol(covariates)
  Ndata = length(BoxCox)
  Nparam = nrow(ParamsGpu)
  
  yx = vclMatrix(cbind(y,
                       matrix(0, n, length(BoxCox)-1),
                       covariates),
                 type=type)      #as.matrix(yx)
  
  coordsGpu<-vclMatrix(spatialmodel$data@coords,type=type)
  
  boxcoxGpu = vclVector(BoxCox, type=type)
  betas = vclMatrix(0, Nparam, Ncov * Ndata, type=type)
  ssqY <- vclMatrix(0, Nparam, Ndata, type=type)
  XVYXVX = vclMatrix(-77, Nparam * Ncov, ncol(yx), type=type)
  ssqBetahat <- vclMatrix(0, Nparam, Ndata, type=type)
  detVar = vclVector(0, Nparam,type=type)
  detReml = vclVector(0, Nparam, type=type)
  jacobian = vclVector(0, Ndata, type=type)   
  ssqYX = vclMatrix(0, ncol(yx) * NparamPerIter, ncol(yx), type=type)
  ssqYXcopy = vclMatrix(0, ncol(yx) * NparamPerIter, ncol(yx), type=type)
  LinvYX = vclMatrix(0, nrow(yx) * NparamPerIter, ncol(yx), type=type)
  QinvSsqYx = vclMatrix(0, NparamPerIter*Ncov, Ndata, type = type)
  cholXVXdiag = vclMatrix(0, NparamPerIter, Ndata, type=type)
  varMat = vclMatrix(0, nrow(yx)*NparamPerIter, nrow(coords), type=type)
  cholDiagMat = vclMatrix(0, NparamPerIter, nrow(coords), type=type)
  
  gpuRandom:::likfitGpu_BackendP(
    yx,        #1
    coordsGpu, #2
    ParamsGpu, #3
    boxcoxGpu,   #4
    betas,  #5
    ssqY,     #6
    XVYXVX,
    ssqBetahat, #8
    detVar,    #9
    detReml,   #10
    jacobian,  #11
    NparamPerIter,  #12
    Nglobal,  #13
    Nlocal,  #14
    NlocalCache,  #15
    verbose = 2,  #16
    ssqYX, 
    ssqYXcopy,  #new
    LinvYX,  #19
    QinvSsqYx, 
    cholXVXdiag, #21
    varMat,        #22     Vbatch
    cholDiagMat)   #new
  
  
  
  # resid^T V^(-1) resid, resid = Y - X betahat = two
  two <- ssqY - ssqBetahat
  minusTwoLogLik <- vclMatrix(0, Nparam, Ndata, type=type)
  temp <- vclMatrix(0, Nparam, Ndata, type=type)
  
  if(form ==2){ #ml
    # = n*log(two/n) + logD + jacobian +n + n*log(2*pi)
    
    matrix_vector_sumBackend(n*log(two/n),
                             detVar,
                             jacobian,  
                             n*log(2*pi),
                             minusTwoLogLik,
                             Nglobal)
    
    
  }else if(form==4){ #mlFixBeta
    # two/variances + n*log(variances) + logD + jacobian + n*log(2*pi)
    
    mat_vec_eledivideBackend(two, variances, temp,  Nglobal)
    
    
    matrix_vector_sumBackend(temp,
                             n*log(variances)+detVar,
                             jacobian,  
                             n*log(2*pi),
                             minusTwoLogLik,
                             Nglobal)
    
    
  }else if(form==5){ #reml
    #first_part <- (n-p)*log(variances) + logD + logP
    #minusTwoLogLik=first_part + two/variances +jacobian + n*log(2*pi)
    
    mat_vec_eledivideBackend(two, variances, temp,  Nglobal)
    
    matrix_vector_sumBackend(temp,
                             detVar+detReml+(n-p)*log(variances),
                             jacobian,  
                             n*log(2*pi),
                             minusTwoLogLik,
                             Nglobal)
    
    
  }else if(form==6){ # remlPro
    # minusTwoLogLik= (n-p)*log(two/(n-p)) + logD + logP + jacobian + n*log(2*pi) + n-p
    
    matrix_vector_sumBackend((n-p)*log(two/(n-p)),
                             detVar+detReml,
                             jacobian,  
                             n*log(2*pi)+n-p,
                             minusTwoLogLik,
                             Nglobal)
    
  }
  
  
}














