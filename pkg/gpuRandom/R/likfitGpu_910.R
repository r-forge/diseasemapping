#' @title Estimate Log-likelihood for Gaussian random fields
#' @useDynLib gpuRandom
#' @export


#before start, we have a spatial model
likfitGpu <- function(spatialmodel,     #data,
                        type = c("double","float"),
                        paramsGpu, #a vclmatrix, consists of all the parameters
                        BoxCox, # an R vector, will always be c(1,0,.....)
                        form = c("loglik", "ml", "mlFixSigma", "mlFixBeta", "reml", "remlPro", "NULL"),
                        NparamPerIter,  # how many sets of params to be evaluated in each loop
                        nBetahats =1,  # how many betahats do you want to calculate for one dataset?
                        minustwotimes=TRUE,
                        Nglobal,
                        Nlocal,
                        NlocalCache,
                        betas=NULL, #a vclmatrix  #only one row batch, but many colbatches
                        verbose=1){

  form = c(loglik=1, ml=2, mlFixSigma=3, mlFixBeta=4, reml=5, remlPro=6, NULL=7)[form]
  if (form==7){
    nBetahats = 0
  }

  covariates = model.matrix(spatialmodel$model$formula, data=spatialmodel$data)
  temp = model.frame(spatialmodel$model$formula, data=spatialmodel$data)
  y = temp[,as.character(attributes(terms(temp))$variables)[2]]
  n = length(y)
  variances <- vclVector(paramsGpu[,3],type=type)     # a vclvector


  # if(BoxCox[1] != 1) {
     BoxCox = c(1, 0, setdiff(BoxCox, c(0,1)))
  # }

  Ncov = ncol(covariates)
  Ndata = length(BoxCox)
  Nparam = nrow(paramsGpu)


  yx = vclMatrix(cbind(y,
                       matrix(0, n, length(BoxCox)-1),
                       covariates),
                 type=type)      #as.matrix(yx)


  coordsGpu<-vclMatrix(spatialmodel$data@coords,type=type)
  boxcoxGpu = vclVector(BoxCox, type=type)
  if(!is.null(betas)){  # do betas computation
    verbose[2]=1
  }else{
    verbose[2]=0
    betasGpu = vclMatrix(0, Ncov, Ndata, type=type)
  }
  ssqY <- vclMatrix(0, Nparam, Ndata, type=type)
  XVYXVX = vclMatrix(0, Nparam * Ncov, ncol(yx), type=type)
  ssqBetahat <- vclMatrix(0, Nparam, Ndata, type=type)
  ssqBeta <- vclMatrix(0, Nparam, Ndata, type=type)
  detVar = vclVector(0, Nparam,type=type)
  detReml = vclVector(0, Nparam, type=type)
  jacobian = vclVector(0, Ndata, type=type)
  ssqYX = vclMatrix(0, ncol(yx) * NparamPerIter, ncol(yx), type=type)
  aTDinvb_beta = vclMatrix(0, Nparam, Ndata, type=type)
  ssqYXcopy = vclMatrix(0, ncol(yx) * NparamPerIter, ncol(yx), type=type)
  LinvYX = vclMatrix(0, nrow(yx) * NparamPerIter, ncol(yx), type=type)
  QinvSsqYx = vclMatrix(0, NparamPerIter*Ncov, Ndata, type = type)
  cholXVXdiag = vclMatrix(0, NparamPerIter, Ncov, type=type)
  varMat = vclMatrix(0, n*NparamPerIter, n, type=type)
  cholDiagMat = vclMatrix(0, NparamPerIter, n, type=type)
  b_beta = vclMatrix(0, NparamPerIter*n, Ndata, type=type)






  gpuRandom:::likfitGpu_BackendP(
    yx,        #1
    coordsGpu, #2
    paramsGpu, #3
    boxcoxGpu,   #4#betasGpu,  #5
    ssqY,     #5#aTDinvb_beta,
    XVYXVX,
    ssqBetahat, #7#ssqBeta,
    detVar,    #11
    detReml,   #12
    jacobian,  #13
    NparamPerIter,  #14
    as.integer(Nglobal),  #12
    as.integer(Nlocal),  #16
    NlocalCache,  #14
    verbose,  #15
    ssqYX, #
    ssqYXcopy,  #new
    LinvYX,  #18
    QinvSsqYx,
    cholXVXdiag, #20
    varMat,        #21     Vbatch
    cholDiagMat) #22
    #b_beta)   #new



  if(form != 7){
  temp <- vclMatrix(0, Nparam, Ndata, type=type)
  # resid^T V^(-1) resid, resid = Y - X betahat = two
  two <- ssqY - ssqBetahat
  minusTwoLogLik <- vclMatrix(0, Nparam, Ndata, type=type)
  }

  if(form ==2){ #ml
    # = n*log(two/n) + logD + jacobian +n + n*log(2*pi)
    matrix_vector_sumBackend(n*log(two/n),
                             detVar,
                             jacobian,
                             n + n*log(2*pi),
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




  ##### calculate betahat #####################
  ###### (X^T V^(-1)X)^(-1) * (X^T V^(-1) y)
  if (nBetahats > 0){
    Betahat <- matrix(0, nrow=nBetahats*Ncov, ncol=Ndata)
    if(nBetahats > 5 | nBetahats > Nparam){
      stop("too many Betahats required")
    }
    for (j in 1:Ndata){
      index <- order(minusTwoLogLik[,j], decreasing = FALSE)[1:nBetahats]  #from small to large
      #print(index)
      for (i in 1:nBetahats){
        a<-c((index[i]-1)*Ncov+1, index[i]*Ncov)
        Betahat[c((i-1)*Ncov+1, i*Ncov),j] <- solve(XVYXVX[a,((Ndata+1):ncol(yx))]) %*% XVYXVX[a,j]
      }
    }
  }





  # if(form ==1 | form ==3){
  #   if(is.null(betas)){
  #     stop("must supply betas for this likelihood")
  #   }else{
  #     one <- ssqY - 2*aTDinvb_beta + ssqBeta
  # 
  #     if(form == 1){
  #       mat_vec_eledivideBackend(one, variances, temp,  Nglobal)
  #       matrix_vector_sumBackend(temp,
  #                                detVar + n*log(variances),
  #                                jacobian,
  #                                n*log(2*pi),
  #                                minusTwoLogLik,
  #                                Nglobal)
  #     }else if(form == 3){
  #       matrix_vector_sumBackend(n*log(one/n),
  #                                detVar,
  #                                jacobian,
  #                                n + n*log(2*pi),
  #                                minusTwoLogLik,
  #                                Nglobal)
  #     }
  #   }
  # }





  if(is.null(betas) & form !=7){
    Theoutput <- list(minusTwoLogLik = minusTwoLogLik,
                      Betahats=Betahat,
                      ssqBetahat = ssqBetahat,
                      ssqY=ssqY,
                      detVar = detVar,   # log |D|
                      detReml = detReml,   # log |P|
                      jacobian = jacobian,
                      XVYXVX=XVYXVX)
  # }else if(form ==1 | form ==3){
  #   Theoutput <- list(minusTwoLogLik = minusTwoLogLik,
  #                     Betahats=Betahat,
  #                     ssqBeta =ssqBeta,
  #                     ssqY=ssqY,
  #                     detVar = detVar,
  #                     detReml = detReml,
  #                     jacobian = jacobian,
  #                     XVYXVX=XVYXVX)
  }else if(form ==7){
    Theoutput <- list(
                      ssqY=ssqY,
                      detVar = detVar,
                      detReml = detReml,
                      jacobian = jacobian,
                      XVYXVX=XVYXVX)
  }


  Theoutput




}















