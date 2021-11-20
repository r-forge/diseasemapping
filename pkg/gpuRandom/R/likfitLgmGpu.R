#' @title Estimate profile Log-likelihood for covariance parameters and lambda
#' @import data.table
#' @useDynLib gpuRandom
#' @export



# betahat and sigmahat 
# profile log-likelihood for each covariance parameters + lambda
likfitLgmGpu <- function(data,
                        formula, 
                        coordinates,
                        paramToEstimate = c("range","shape","nugget","anisoAngleDegrees", "anisoRatio"), 
                        #variance and regression parameters are always estimated
                        params=NULL, # CPU matrix for now, users need to provide proper parameters given their specific need
                        BoxCox = c(1, 0, 0.5),
                        type = c("double","float"),
                        reml=FALSE, 
                        minustwotimes=TRUE,
                        NparamPerIter,
                        Nglobal,
                        Nlocal,
                        NlocalCache,
                        verbose=FALSE){
  
       # will always be c(1,0,.......)
       BoxCox = c(1, 0, setdiff(BoxCox, c(0,1)))
       
       if(length(paramToEstimate)>2)
       {
         warning("change param matrix for another profile log-likelihood\n")
       }
       
       
       
       # get rid of NAs in data
       data = data.frame(data)
       theNA = apply(  data[,all.vars(formula),drop=FALSE],
                        1, 
                       function(qq) any(is.na(qq))
                    )
       noNA = !theNA
       
       
       ############## the X matrix #################################
       covariates = model.matrix(formula, data[noNA,])
       observations = all.vars(formula)[1]
       observations = data.matrix(data[noNA, observations, drop=FALSE])
  
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
       
       #betas <- matrix(0,nrow=Ncov, ncol=Ndata)
       #betasGpu = vclMatrix(betas, type=type)
       ssqY <- vclMatrix(0, Nparam, Ndata, type=type)
       XVYXVX = vclMatrix(0, Nparam * Ncov, ncol(yx), type=type)
       ssqBetahat <- vclMatrix(0, Nparam, Ndata, type=type)
       ssqBeta <- vclMatrix(0, Nparam, Ndata, type=type)
       detVar = vclVector(0, Nparam,type=type)
       detReml = vclVector(0, Nparam, type=type)
       jacobian = vclVector(0, Ndata, type=type)   
       ssqYX = vclMatrix(0, ncol(yx) * NparamPerIter, ncol(yx), type=type)
       #aTDinvb_beta = vclMatrix(0, Nparam, Ndata, type=type)
       ssqYXcopy = vclMatrix(0, ncol(yx) * NparamPerIter, ncol(yx), type=type)
       LinvYX = vclMatrix(0, nrow(yx) * NparamPerIter, ncol(yx), type=type)
       QinvSsqYx = vclMatrix(0, NparamPerIter*Ncov, Ndata, type = type)
       cholXVXdiag = vclMatrix(0, NparamPerIter, Ncov, type=type)
       varMat = vclMatrix(0, Nobs*NparamPerIter, Nobs, type=type)
       cholDiagMat = vclMatrix(0, NparamPerIter, Nobs, type=type)
       #b_beta = vclMatrix(0, NparamPerIter*Nobs, Ndata, type=type)
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
         cholDiagMat)
       
       # any(is.na(as.vector(detVar)))
       # any(is.na(as.matrix(varMat)))
       # any(is.na(as.matrix(ssqYX)))
       # any(is.na(as.matrix(ssqYXcopy)))
       # any(is.na(as.matrix(ssqBetahat)))
       # any(is.na(as.vector(detReml)))
       # any(is.na(as.matrix(cholXVXdiag)))
       # 
       # 
       # any(is.na(as.matrix(Nobs*log(two/Nobs))))
       
       # resid^T V^(-1) resid, resid = Y - X betahat = ssqResidual
       ssqResidual <- ssqY - ssqBetahat
       #any(is.na(as.matrix(ssqResidual)))
       
       
    
       
       
       
       if(reml== FALSE){ # ml
       gpuRandom:::matrix_vector_sumBackend(Nobs*log(ssqResidual/Nobs),
                                detVar,
                                jacobian,  
                                Nobs + Nobs*log(2*pi),
                                minusTwoLogLik,
                                Nglobal)
       }else if(reml==TRUE){ # remlpro
         # minusTwoLogLik= (n-p)*log(two/(n-p)) + logD + logP + jacobian + n*log(2*pi) + n-p    
         gpuRandom:::matrix_vector_sumBackend((Nobs-Ncov)*log(ssqResidual/(Nobs-Ncov)),
                                  detVar+detReml,
                                  jacobian,  
                                  Nobs*log(2*pi)+Nobs-Ncov,
                                  minusTwoLogLik,
                                  Nglobal)
         
       }

       
       LogLikcpu <- as.matrix(-0.5*minusTwoLogLik)
       
       
       
       ##############output matrix####################
       Table <- matrix(NA, nrow=length(paramToEstimate)+Ncov+2, ncol=3)
       rownames(Table) <-  c(paste(c('betahat'), seq_len(Ncov),sep = ''), "sdSpatial", paramToEstimate,  "BoxCox")
       colnames(Table) <-  c("estimate", "Lower95ci", "Upper95ci") #"95Lowerci", "95Upperci")
       
       
  
       
       
      ###############profile for covariance parameters#####################
       if(is.element('range',paramToEstimate)){
         # get profile log-lik for range values
         result = data.table::as.data.table(cbind(LogLikcpu, params0[,"range"]))
         colnames(result) <- c(paste(c('boxcox'),BoxCox,sep = ''), "range")
         profileLogLik <- result[, .(profile=max(.SD)), by=range]
 
         #maximum <- max(profileLogLik$profile)
         #MLE <- profileLogLik$range[which.max(profileLogLik$profile)]
         # plot(profileLogLik$range, profileLogLik$profile)
         f1<-approxfun(profileLogLik$range, profileLogLik$profile)  
         # curve(f1, add=TRUE)
         optimization <- optimize(f1, interval = c(min(profileLogLik$range), max(profileLogLik$range)),maximum = TRUE)
         MLE <- optimization$maximum 
         maximum <- optimization$objective
         
         leftOfMax = profileLogLik$range < profileLogLik$range[which.max(profileLogLik$profile)]
         if(length(which(leftOfMax)) <=2 | length(which(!leftOfMax)) <=2){
           ci95 = c(NA,NA)
           print("Not enough data for CI calculation")
         }else{
         afLeft <- approxfun(profileLogLik$profile[leftOfMax], profileLogLik$range[leftOfMax])   
         # plot(profileLogLik$profile[leftOfMax], profileLogLik$range[leftOfMax])# curve(afLeft, add=TRUE)
         afRight <- approxfun(profileLogLik$profile[!leftOfMax], profileLogLik$range[!leftOfMax])   # plot(profileLogLik$range, profileLogLik$profile)# curve(af, add=TRUE)
         # breaks = maximum - qchisq(0.95,  df = 1)/2
         # ci = c(afLeft(breaks), afRight(breaks))
         # if (any(is.na(ci))){
         #   warning("'range' scope too small for 95% ci")
         # }
         # 
         # ci<-uniroot(function(x) {af(x)-breaks}, interval=c(1,10)
         #         )$root
         breaks95 = maximum - qchisq(0.95,  df = 1)/2
         ci95= c(afLeft(breaks95), afRight(breaks95))
         if (any(is.na(ci95))){
              warning("'range' scope too small for 95% ci")
           }
         }
         Table["range",] <- c(MLE, ci95)
        }
       
       
       
       
       
       if(is.element('shape',paramToEstimate)){
         result = as.data.table(cbind(LogLikcpu, params0[,"shape"]))
         colnames(result) <- c(paste(c('boxcox'),BoxCox,sep = ''), "shape")
         profileLogLik <- result[, .(profile=max(.SD)), by=shape]

         maximum <- max(profileLogLik$profile)
         MLE <- profileLogLik$shape[which.max(profileLogLik$profile)]
         
         leftOfMax = profileLogLik$shape < profileLogLik$shape[which.max(profileLogLik$profile)]
         afLeft <- approxfun(profileLogLik$profile[leftOfMax], profileLogLik$shape[leftOfMax])  
         afRight <- approxfun(profileLogLik$profile[!leftOfMax], profileLogLik$shape[!leftOfMax])   
         # breaks95 = maximum - qchisq(0.95,  df = 1)/2
         # ci95= c(afLeft(breaks95), afRight(breaks95))

         #plot(profileLogLik$shape, profileLogLik$profile)
         
         breaks95 = maximum - qchisq(0.95,  df = 1)/2
         ci95= c(afLeft(breaks95), afRight(breaks95))
         if (any(is.na(ci95)))
           warning("'shape' scope too small for 95% ci")
         
         
         Table["shape",] <- c(MLE,  ci95)
       }       
  
  
       if(is.element('nugget',paramToEstimate)){
         result = as.data.table(cbind(LogLikcpu, params0[,"nugget"]))
         colnames(result) <- c(paste(c('boxcox'),boxcox,sep = ''), "nugget")

         profileLogLik <-result[, .(profile=max(.SD)), by=nugget]
         maximum <- max(profileLogLik$profile)
         MLE <- profileLogLik$nugget[which.max(profileLogLik$profile)]
         
         leftOfMax = profileLogLik$nugget < profileLogLik$nugget[which.max(profileLogLik$profile)]
         afLeft <- approxfun(profileLogLik$profile[leftOfMax], profileLogLik$nugget[leftOfMax])   # plot(profileLogLik$range, profileLogLik$profile)# curve(af, add=TRUE)
         afRight <- approxfun(profileLogLik$profile[!leftOfMax], profileLogLik$nugget[!leftOfMax])   # plot(profileLogLik$range, profileLogLik$profile)# curve(af, add=TRUE)
         breaks95 = maximum - qchisq(0.95,  df = 1)/2
         ci95= c(afLeft(breaks95), afRight(breaks95))
         if (any(is.na(ci95)))
           warning("'nugget' scope too small for 95% ci")
         
         Table["nugget",] <- c(MLE,  ci95)
       }
  
       
       
       
       if(is.element('anisoRatio',paramToEstimate)){
         result = as.data.table(cbind(LogLikcpu, params0[,"anisoRatio"]))
         colnames(result) <- c(paste(c('boxcox'),boxcox,sep = ''), "anisoRatio")
         profileLogLik <- result[, .(profile=max(.SD)), by=anisoRatio]
         maximum <- max(profileLogLik$profile)
         MLE <- profileLogLik$anisoRatio[which.max(profileLogLik$profile)]
         
         leftOfMax = profileLogLik$anisoRatio < profileLogLik$anisoRatio[which.max(profileLogLik$profile)]
         afLeft <- approxfun(profileLogLik$profile[leftOfMax], profileLogLik$anisoRatio[leftOfMax])   # plot(profileLogLik$range, profileLogLik$profile)# curve(af, add=TRUE)
         afRight <- approxfun(profileLogLik$profile[!leftOfMax], profileLogLik$anisoRatio[!leftOfMax])   # plot(profileLogLik$range, profileLogLik$profile)# curve(af, add=TRUE)
         
         breaks95 = maximum - qchisq(0.95,  df = 1)/2
         ci95= c(afLeft(breaks95), afRight(breaks95))
         if (any(is.na(ci95)))
           warning("'anisoRatio' scope too small for 95% ci")
         
         Table["anisoRatio",] <- c(MLE,  ci95)
       }
       
       
       
       if(is.element('anisoAngleDegrees',paramToEstimate)){
         result = as.data.table(cbind(LogLikcpu, params0[,"anisoAngleDegrees"]))
         colnames(result) <- c(paste(c('boxcox'),boxcox,sep = ''), "anisoAngleDegrees")
         profileLogLik <- result[, .(profile=max(.SD)), by=anisoAngleDegrees]
         
         maximum <- max(profileLogLik$profile)
         MLE <- profileLogLik$anisoAngleDegrees[which.max(profileLogLik$profile)]
         
         leftOfMax = profileLogLik$anisoAngleDegrees < profileLogLik$anisoAngleDegrees[which.max(profileLogLik$profile)]
         afLeft <- approxfun(profileLogLik$profile[leftOfMax], profileLogLik$anisoAngleDegrees[leftOfMax])   # plot(profileLogLik$range, profileLogLik$profile)# curve(af, add=TRUE)
         afRight <- approxfun(profileLogLik$profile[!leftOfMax], profileLogLik$anisoAngleDegrees[!leftOfMax])   # plot(profileLogLik$range, profileLogLik$profile)# curve(af, add=TRUE)
         
         breaks95 = maximum - qchisq(0.95,  df = 1)/2
         ci95= c(afLeft(breaks95), afRight(breaks95))
         if (any(is.na(ci95)))
           warning("'anisoAngleDegrees' scope too small for 95% ci")
         
         Table["anisoAngleDegrees",] <- c(MLE,  ci95)
       }
       
  
       
       
       
       ###############lambda hat#####################
    
       #likForBoxCox = apply( LogLikcpu, 2, max )
       index <- which(LogLikcpu == max(LogLikcpu, na.rm = TRUE), arr.ind = TRUE)
       
       Table["BoxCox",1] <- BoxCox[index[2]]
       
       #any(is.na(as.matrix(minusTwoLogLik)))
       
       ###############betahat#####################
       Betahat <- matrix(0, nrow=Ncov, ncol=Ndata)
       a<-c((index[1]-1)*Ncov+1, index[1]*Ncov)
       Betahat <- solve(XVYXVX[a,((Ndata+1):ncol(yx))]) %*% XVYXVX[a,index[2]]
       
       Table[paste(c('betahat'),seq_len(Ncov),sep = ''),1] <- Betahat
       
       
       #################sigma hat#########################
       if(reml==FALSE)  {
         Table["sdSpatial",1] <- sqrt(ssqResidual[index[1],index[2]]/Nobs)
       }
       else if(reml==TRUE) {         
         Table["sdSpatial",1] <- sqrt(ssqResidual[index[1],index[2]]/(Nobs - Ncov))
       }
       


       Output <- list(LogLik=LogLikcpu,
                      minusTwoLogLik = minusTwoLogLik,
                      table = Table,
                      ssqY=ssqY,     
                      ssqBetahat = ssqBetahat,
                      detVar = detVar,   
                      detReml = detReml,   
                      jacobian = jacobian,
                      XVYXVX=XVYXVX)
       
       Output

       
       
       
       
       
       
       
       
       
       
       
}

