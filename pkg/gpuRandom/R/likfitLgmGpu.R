#' @title Estimate Log-likelihood for Gaussian random fields
#' @useDynLib gpuRandom
#' @export












library(data.table)

# betahat and sigmahat 
# profile log-likelihood for each covariance parameters + lambda
likfitLgmGpu <- function(data,
                        formula, 
                        coordinates=data,
                        paramToEstimate = c("range","shape","nugget","anisoAngleDegrees", "anisoRatio"), 
                        #variance and regression parameters are always estimated
                        params=NULL, # CPU matrix for now, users need to provide proper parameters given their specific need
                        BoxCox = c(1, 0, 0.5),
                        type = c("double","float"),
                        reml=TRUE, 
                        minustwotimes=TRUE,
                        NparamPerIter,
                        upper=NULL,
                        lower=NULL, 
                        Nglobal,
                        Nlocal,
                        NlocalCache,
                        verbose=FALSE){
  
       # will always be c(1,0,.......)
       BoxCox = c(1, 0, setdiff(BoxCox, c(0,1)))
       
       # get rid of NAs in data
       data = data.frame(data)
       theNA = apply(  data[,all.vars(formula),drop=FALSE],
                        1, 
                       function(qq) any(is.na(qq))
                    )
       noNA = !theNA
       
       covariates = model.matrix(formula, data[noNA,])
       observations = all.vars(formula)[1]
       observations = data[noNA, observations, drop=FALSE]
  
       Nobs = nrow(data)
       Ncov = ncol(covariates)
       Ndata = length(BoxCox)
       Nparam = nrow(params)

       # whole data set including columns for transformed y's
       yx = vclMatrix(cbind(observations, matrix(0, Nobs, Ndata-1), covariates), 
                      type=type)
       

       
       # coordinates
       coordsGpu = vclMatrix(coordinates, type=type)  
       # box-cox
       boxcoxGpu = vclVector(BoxCox, type=type)
       
       # prepare 
       params0 = geostatsp::fillParam(params)
       paramsGpu = vclMatrix(cbind(params0, matrix(0, nrow(params0), 22-ncol(params0))),type=type)
       gpuR::colnames(paramsGpu) = colnames(params0)
       
       
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
       varMat = vclMatrix(0, Nobs*NparamPerIter, Nobs, type=type)
       cholDiagMat = vclMatrix(0, NparamPerIter, Nobs, type=type)
       b_beta = vclMatrix(0, NparamPerIter*Nobs, Ndata, type=type)
       minusTwoLogLik <- vclMatrix(0, Nparam, Ndata, type=type)
  
       # paramDefaults = c(
       #   nugget=0,
       #   anisoRatio=1, 
       #   anisoAngleRadians=0,
       #   shape=1.5, boxcox=1,
       #   range=maxDist/10
       # )
       
       gpuRandom:::likfitGpu_BackendP(
         yx,        #1
         coordsGpu, #2
         paramsGpu, #3
         boxcoxGpu,   #4
         betasGpu,  #5
         ssqY,     #6
         aTDinvb_beta,
         XVYXVX,
         ssqBetahat, #9
         ssqBeta,
         detVar,    #11
         detReml,   #12
         jacobian,  #13
         NparamPerIter,  #14
         as.integer(Nglobal),  #15
         as.integer(Nlocal),  #16
         NlocalCache,  #17
         verbose,  #18
         ssqYX, #
         ssqYXcopy,  #new
         LinvYX,  #19
         QinvSsqYx, 
         cholXVXdiag, #22
         varMat,        #23     Vbatch
         cholDiagMat,
         b_beta)   #new 245
       

       
       # resid^T V^(-1) resid, resid = Y - X betahat = two
       two <- ssqY - ssqBetahat
       
 
       
       gpuRandom:::matrix_vector_sumBackend(n*log(two/n),
                                detVar,
                                jacobian,  
                                n + n*log(2*pi),
                                minusTwoLogLik,
                                Nglobal)
       
       
        Betahat <- matrix(0, nrow=Ncov, ncol=Ndata)
        LogLikcpu <- as.matrix(-0.5*minusTwoLogLik)
        index <- which(LogLikcpu == max(LogLikcpu, na.rm = TRUE), arr.ind = TRUE)

        a<-c((index[1]-1)*Ncov+1, index[i]*Ncov)
        Betahat[c((i-1)*Ncov+1, i*Ncov),index[2]] <- solve(XVYXVX[a,((Ndata+1):ncol(yx))]) %*% XVYXVX[a,index[2]]
         
 
       
       output <- matrix(NA, nrow=7, ncol=4)
       rownames(output) <-  c("range","shape","nugget","anisoRatio", "anisoAngleDegrees", "variance", "betas")
       colnames(output) <-  c("estimate", "LogLik", "95Lowerci", "95Upperci")
       
       
       
       
       if(is.element('range',paramToEstimate)){
         result = as.data.table(cbind(LogLikcpu, params0[,"range"]))
         colnames(result) <- c(paste(c('boxcox'),boxcox,sep = ''), "range")
         profileLogLik <- result[, .(profile=max(.SD)), by=range]
         maximum <- max(profileLogLik$profile)
         MLE_range <- profileLogLik$range[which.max(profileLogLik$profile)]
         
         af <- approxfun(profileLogLik$range, profileLogLik$profile)   # plot(profileLogLik$range, profileLogLik$profile)# curve(af, add=TRUE)
         breaks = maximum - qchisq(0.95,  df = 1)/2
         ci<-uniroot(function(x) {af(x)-breaks}, interval=c(1,10)
                 )$root

         output["range",] <- c(MLE_range, maximum, ci, NA)
        }
       
       if(is.element('shape',paramToEstimate)){
         result = as.data.table(cbind(as.matrix(-0.5*minusTwoLogLik), params0[,"shape"]))
         colnames(result) <- c(paste(c('boxcox'),boxcox,sep = ''), "shape")
         result[, .(profile=max(.SD)), by=shape]
         
       }       
  
  
       if(is.element('nugget',paramToEstimate)){
         result = as.data.table(cbind(as.matrix(-0.5*minusTwoLogLik), params0[,"nugget"]))
         colnames(result) <- c(paste(c('boxcox'),boxcox,sep = ''), "nugget")
         result[, .(profile=max(.SD)), by=nugget]
         
       }
  
       if(is.element('anisoRatio',paramToEstimate)){
         result = as.data.table(cbind(as.matrix(-0.5*minusTwoLogLik), params0[,"anisoRatio"]))
         colnames(result) <- c(paste(c('boxcox'),boxcox,sep = ''), "anisoRatio")
         result[, .(profile=max(.SD)), by=anisoRatio]
         
       }
       
       
       
       if(is.element('anisoAngleRadians',paramToEstimate)){
         result = as.data.table(cbind(as.matrix(-0.5*minusTwoLogLik), params0[,"anisoAngleRadians"]))
         colnames(result) <- c(paste(c('boxcox'),boxcox,sep = ''), "anisoAngleRadians")
         result[, .(profile=max(.SD)), by=anisoAngleRadians]
         
       }
       
       
       if(is.element('anisoAngleDegrees',paramToEstimate)){
         result = as.data.table(cbind(as.matrix(-0.5*minusTwoLogLik), params0[,"anisoAngleDegrees"]))
         colnames(result) <- c(paste(c('boxcox'),boxcox,sep = ''), "anisoAngleDegrees")
         result[, .(profile=max(.SD)), by=anisoAngleDegrees]
         
       }
       
       
       
       
       
       
       
       
       
       
       
       
       
}

