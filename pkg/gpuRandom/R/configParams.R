#' @title set loglikelihood locations
#' @useDynLib gpuRandom
#' @export


    getHessian <- function(Model,
                           Mle = NULL){
  
      delta = 0.01  
      kappa = 60
      # order the params
      Theorder <- c('range', 'shape', 'nugget', 'anisoRatio', 'anisoAngleRadians')
      parToLog <- c("range", 'shape', "nugget")
      if(is.null(Mle)){
        Mle <- Model$optim$options[,'opt']
      }
      Mle <- Mle[order(match(names(Mle), Theorder))]
      Mle <- Mle[!names(Mle) %in% c('boxcox')]
      
      ## check mle nugget    
      if(Mle['nugget'] < delta){
        Mle['nugget'] = delta 
        parToLog <- parToLog[!parToLog %in% 'nugget']
      }  
      
      whichLogged = which(names(Mle) %in% parToLog)
      whichAniso = which(names(Mle) %in% c('anisoRatio', 'anisoAngleRadians'))
      
      if('anisoRatio' %in% names(Mle)){
        if(!('anisoAngleRadians' %in% names(Mle))){
          stop('anisoRatio and anisoAngleRadians must be together')
        }
        
        if(Mle['anisoRatio'] <= 1){
          gamma3 <-  unname(sqrt(1/Mle['anisoRatio']-1) * cos(2*Mle['anisoAngleRadians']))
          gamma4 <-  unname(sqrt(1/Mle['anisoRatio']-1) * sin(2*Mle['anisoAngleRadians']))
        }else{
          gamma3 <-  unname(sqrt(Mle['anisoRatio']-1) * cos(2*Mle['anisoAngleRadians']))
          gamma4 <-  unname(sqrt(Mle['anisoRatio']-1) * sin(2*Mle['anisoAngleRadians']))
        }
        aniso <- c(gamma3 = gamma3, gamma4 = gamma4)
        if(whichLogged[2]-whichLogged[1]>1){
          MleGamma = c(log(Mle[whichLogged[1]]), Mle[-c(whichLogged, whichAniso)],log(Mle[whichLogged[2]]), aniso) 
        }else{
          MleGamma = c(log(Mle[whichLogged]), Mle[-c(whichLogged, whichAniso)], aniso)
        }
      }else{
        MleGamma = c(log(Mle[whichLogged]),  Mle[-whichLogged])
      }
      names(MleGamma)[whichLogged] = paste("log(", names(Mle)[whichLogged], ")",sep="")
      
  ## center approximation
  a1 <- c(1,1,-1,-1)
  a2 <- c(1,-1,1,-1)
  a0 <- rep(0,4)
  if(length(Mle)==1){
    derivGridDf <- matrix(c(-1,0,1), nrow=3, ncol=1)
    
  }else if(length(Mle)==2){
    diagonals <- rbind(rep(0,4),
                       c(1,0,0,0),
                       c(-1,0,0,0),
                       c(0, 1,0,0),
                       c(0,-1,0,0))
    
    derivGridDf <- rbind(
      diagonals,
      cbind(a1, a2, a0, a0))[,1:2]
    
  }else if(length(Mle)==3){
    diagonals <- rbind(rep(0,3),
                       c(1,0,0),
                       c(-1,0,0),
                       c(0,1,0),
                       c(0,-1,0),
                       c(0,0,1),
                       c(0,0,-1))
    derivGridDf <- rbind(
      diagonals,
      cbind(a1, a2, a0),
      cbind(a1, a0, a2),
      cbind(a0, a1, a2))
    
    
  }else if(length(Mle)==4){
    diagonals <- rbind(rep(0,4),
                       c(1,0,0,0),
                       c(-1,0,0,0),
                       c(0,1,0,0),
                       c(0,-1,0,0),
                       c(0,0,1,0),
                       c(0,0,-1,0),
                       c(0,0,0,1),
                       c(0,0,0,-1))
    derivGridDf <- rbind(
      diagonals,
      cbind(a1, a2, a0, a0),
      cbind(a1, a0, a2, a0),
      cbind(a1, a0, a0, a2),
      cbind(a0, a1, a2, a0),
      cbind(a0, a1, a0, a2),
      cbind(a0, a0, a1, a2))
  }else if(length(Mle)==5){
    diagonals <- rbind(rep(0,5),
                       c(1,0,0,0,0),
                       c(-1,0,0,0,0),
                       c(0,1,0,0,0),
                       c(0,-1,0,0,0),
                       c(0,0,1,0,0),
                       c(0,0,-1,0,0),
                       c(0,0,0,1,0),
                       c(0,0,0,-1,0),
                       c(0,0,0,0,1),
                       c(0,0,0,0,-1))
    
    derivGridDf <- rbind(
      diagonals,
      cbind(a1, a2, a0, a0, a0),
      cbind(a1, a0, a2, a0, a0),
      cbind(a1, a0, a0, a2, a0),
      cbind(a1, a0, a0, a0, a2),
      cbind(a0, a1, a2, a0, a0),
      cbind(a0, a1, a0, a2, a0),
      cbind(a0, a1, a0, a0, a2),
      cbind(a0, a0, a1, a2, a0),
      cbind(a0, a0, a1, a0, a2),
      cbind(a0, a0, a0, a1, a2))
  }
  
  deltas = rep(0.01, length(Mle))
  # names(deltas) = names(MleGamma)
  # deltas[c('log(shape)','log(nugget)')] = 1e-3
  if(length(Mle)<=1){
    ParamsetGamma <- matrix(MleGamma, nrow=nrow(derivGridDf), ncol=length(Mle), byrow=TRUE, dimnames = list(NULL, names(MleGamma))) + derivGridDf*deltas
  }else{
    ParamsetGamma <- matrix(MleGamma, nrow=nrow(derivGridDf), ncol=length(Mle), byrow=TRUE, dimnames = list(NULL, names(MleGamma))) + 
      derivGridDf %*% diag(deltas)
  }
  
  if(!('anisoRatio' %in% names(Mle))){
    Paramset <- cbind(exp(ParamsetGamma[,paste("log(", names(Mle)[whichLogged], ")",sep="")]), ParamsetGamma[,-whichLogged])   
  }else{
    #    if(Mle['anisoRatio'] <= 1) {
    # Paramset <- cbind(exp(ParamsetGamma[,paste("log(", parToLog, ")",sep="")]), 
    #                         1/(sqrt(ParamsetGamma[,'gamma3']^2 + ParamsetGamma[,'gamma4']^2) + 1), atan(ParamsetGamma[,'gamma4']/ParamsetGamma[,'gamma3'])/2)        
    #    }else{
    temp <- as.data.frame(ParamsetGamma[,'gamma3'] + 1i * ParamsetGamma[,'gamma4'])
    naturalspace <- cbind(Mod(temp[,1])^2 + 1, Arg(temp[,1])/2)
    if(whichLogged[2]-whichLogged[1]>1){
      Paramset <- cbind(exp(ParamsetGamma[, whichLogged[1]]), ParamsetGamma[,-c(whichLogged, whichAniso)], exp(ParamsetGamma[, whichLogged[2]]), naturalspace)
    }else{
      Paramset <- cbind(exp(ParamsetGamma[, whichLogged]), ParamsetGamma[,-c(whichLogged, whichAniso)], naturalspace)
    }
  } 
  
  colnames(Paramset) <- names(Mle)  
  toAdd = setdiff(c('range','shape','nugget','anisoRatio', 'anisoAngleRadians'), names(Mle))
  otherParams = matrix(Model$opt$mle[toAdd], nrow=nrow(Paramset), ncol = length(toAdd),
                       dimnames = list(rownames(Paramset), toAdd), byrow=TRUE)
  Params1 <- cbind(Paramset, otherParams)
  
  
  
  result1<-gpuRandom::getProfLogL(data= Model$data,
                                 formula=Model$model$formula,
                                 coordinates=Model$data@coords,
                                 params=Params1,
                                 boxcox = Model$parameters['boxcox'],
                                 type = "double",
                                 NparamPerIter=100,
                                 gpuElementsOnly=FALSE,
                                 reml=FALSE,
                                 Nglobal=c(128,64),
                                 Nlocal=c(16,16),
                                 NlocalCache=2000)
  
  
  #result$paramsRenew
  A <- result1$LogLik[-1, 2]
  Origin <- result1$LogLik[1, 2]
  #plot(result$LogLik)
  
  HessianMat <- matrix(0, nrow=length(Mle), ncol=length(Mle))
  rownames(HessianMat) <- names(MleGamma)
  colnames(HessianMat) <- names(MleGamma)
  
  
  if(length(Mle)==1){
    HessianMat[1,1] <- (result$LogLik[1, 2] - 2*result$LogLik[2, 2] + result$LogLik[3, 2])/(deltas[1]^2)
  }else if(length(Mle)==2){
    HessianMat[1,1] <- (A[1] + A[2] -2*Origin)/(deltas[1]^2)
    HessianMat[2,2] <- (A[3] + A[4] -2*Origin)/(deltas[2]^2)
    HessianMat[1,2] <- (A[5] - A[6] - A[7] + A[8])/(4*deltas[1]*deltas[2])
    HessianMat[2,1] <- HessianMat[1,2]
  }else if(length(Mle)==3){
    HessianMat[1,1] <- (A[1] + A[2] -2*Origin)/(deltas[1]^2)
    HessianMat[2,2] <- (A[3] + A[4] -2*Origin)/(deltas[2]^2)
    HessianMat[3,3] <- (A[5] + A[6] -2*Origin)/(deltas[3]^2)
    HessianMat[1,2] <- (A[7] - A[8] - A[9] + A[10])/(4*deltas[1]*deltas[2])
    HessianMat[1,3] <- (A[11] - A[12] - A[13] + A[14])/(4*deltas[1]*deltas[3])
    HessianMat[2,3] <- (A[15] - A[16] - A[17] + A[18])/(4*deltas[2]*deltas[3])
    HessianMat[2,1] <- HessianMat[1,2]
    HessianMat[3,1] <- HessianMat[1,3]
    HessianMat[3,2] <- HessianMat[2,3]
  }else if(length(Mle)==4){
    HessianMat[1,1] <- (A[1] + A[2] -2*Origin)/(deltas[1]^2)
    HessianMat[2,2] <- (A[3] + A[4] -2*Origin)/(deltas[2]^2)
    HessianMat[3,3] <- (A[5] + A[6] -2*Origin)/(deltas[3]^2)
    HessianMat[4,4] <- (A[7] + A[8] -2*Origin)/(deltas[4]^2)
    
    HessianMat[1,2] <- (A[9] - A[10] - A[11] + A[12])/(4*deltas[1]*deltas[2])
    HessianMat[1,3] <- (A[13] - A[14] - A[15] + A[16])/(4*deltas[1]*deltas[3])
    HessianMat[1,4] <- (A[17] - A[18] - A[19] + A[20])/(4*deltas[1]*deltas[4])
    HessianMat[2,3] <- (A[21] - A[22] - A[23] + A[24])/(4*deltas[2]*deltas[3])
    HessianMat[2,4] <- (A[25] - A[26] - A[27] + A[28])/(4*deltas[2]*deltas[4])
    HessianMat[3,4] <- (A[29] - A[30] - A[31] + A[32])/(4*deltas[3]*deltas[4])
    
  }else if(length(Mle)==5){
    HessianMat[1,1] <- (A[1] + A[2] -2*Origin)/(deltas[1]^2)
    HessianMat[2,2] <- (A[3] + A[4] -2*Origin)/(deltas[2]^2)
    HessianMat[3,3] <- (A[5] + A[6] -2*Origin)/(deltas[3]^2)
    HessianMat[4,4] <- (A[7] + A[8] -2*Origin)/(deltas[4]^2)
    HessianMat[5,5] <- (A[9] + A[10] -2*Origin)/(deltas[5]^2)
    
    HessianMat[1,2] <- (A[11] - A[12] - A[13] + A[14])/(4*deltas[1]*deltas[2])
    HessianMat[1,3] <- (A[15] - A[16] - A[17] + A[18])/(4*deltas[1]*deltas[3])
    HessianMat[1,4] <- (A[19] - A[20] - A[21] + A[22])/(4*deltas[1]*deltas[4])
    HessianMat[1,5] <- (A[23] - A[24] - A[25] + A[26])/(4*deltas[1]*deltas[5])
    HessianMat[2,3] <- (A[27] - A[28] - A[29] + A[30])/(4*deltas[2]*deltas[3])
    HessianMat[2,4] <- (A[31] - A[32] - A[33] + A[34])/(4*deltas[2]*deltas[4])
    HessianMat[2,5] <- (A[35] - A[36] - A[37] + A[38])/(4*deltas[2]*deltas[5])
    HessianMat[3,4] <- (A[39] - A[40] - A[41] + A[42])/(4*deltas[3]*deltas[4])
    HessianMat[3,5] <- (A[43] - A[44] - A[45] + A[46])/(4*deltas[3]*deltas[5])
    HessianMat[4,5] <- (A[47] - A[48] - A[49] + A[50])/(4*deltas[4]*deltas[5])
    
    HessianMat[5,1] <- HessianMat[1,5]
    HessianMat[5,2] <- HessianMat[2,5]
    HessianMat[5,3] <- HessianMat[3,5]
    HessianMat[5,4] <- HessianMat[4,5]
  }
  
  if(length(Mle)==4 | length(Mle)==5){
    HessianMat[2,1] <- HessianMat[1,2]
    HessianMat[3,1] <- HessianMat[1,3]
    HessianMat[3,2] <- HessianMat[2,3]
    HessianMat[4,1] <- HessianMat[1,4]
    HessianMat[4,2] <- HessianMat[2,4]
    HessianMat[4,3] <- HessianMat[3,4]
  }

  
  
  ## get the First derivative
  # check if shape and nugget's mle is close to 0
  # whichLogged = which(names(Mle) %in% parToLog)
  # whichAniso = which(names(Mle) %in% c('anisoRatio', 'anisoAngleRadians'))
  # 
  # if('anisoRatio' %in% names(Mle)){
  #   if(!('anisoAngleRadians' %in% names(Mle))){
  #     stop('anisoRatio and anisoAngleRadians must be together')
  #   }
  #   
  #   if(Mle['anisoRatio'] <= 1){
  #     gamma3 <-  unname(sqrt(1/Mle['anisoRatio']-1) * cos(2*Mle['anisoAngleRadians']))
  #     gamma4 <-  unname(sqrt(1/Mle['anisoRatio']-1) * sin(2*Mle['anisoAngleRadians']))
  #   }else{
  #     gamma3 <-  unname(sqrt(Mle['anisoRatio']-1) * cos(2*Mle['anisoAngleRadians']))
  #     gamma4 <-  unname(sqrt(Mle['anisoRatio']-1) * sin(2*Mle['anisoAngleRadians']))
  #   }
  #   aniso <- c(gamma3 = gamma3, gamma4 = gamma4)
  #   if(whichLogged[2]-whichLogged[1]>1){
  #     MleGamma = c(log(Mle[whichLogged[1]]), Mle[-c(whichLogged, whichAniso)],log(Mle[whichLogged[2]]), aniso) 
  #   }else{
  #     MleGamma = c(log(Mle[whichLogged]), Mle[-c(whichLogged, whichAniso)], aniso)
  #   }
  # }else{
  #   MleGamma = c(log(Mle[whichLogged]),  Mle[-whichLogged])
  # }
  # names(MleGamma)[whichLogged] = paste("log(", names(Mle)[whichLogged], ")",sep="")
  
  
  # frst detivatives
  derivGridDf1 <- rbind(c(1,0,0,0,0),
                        c(-1, 0,0,0,0),
                        c(0, 1, 0,0,0),
                        c(0,-1, 0,0,0),
                        c(0, 0, 1,0,0),
                        c(0, 0,-1,0,0),
                        c(0, 0, 0,1,0),
                        c(0, 0, 0,-1,0),
                        c(0, 0, 0, 0,1),
                        c(0, 0, 0, 0,-1))
  
  deltas = rep(delta, length(Mle))
  # names(deltas) = names(MleGamma)
  # deltas['gamma4'] = 0.02
  ParamsetGamma <- matrix(MleGamma, nrow=nrow(derivGridDf1), ncol=length(Mle), byrow=TRUE, dimnames = list(NULL, names(MleGamma))) + 
    derivGridDf1 %*% diag(deltas)
  
  if(!('anisoRatio' %in% names(Mle))){
    Paramset <- cbind(exp(ParamsetGamma[,paste("log(", names(Mle)[whichLogged], ")",sep="")]), ParamsetGamma[,-whichLogged])   
  }else{
    temp <- as.data.frame(ParamsetGamma[,'gamma3'] + 1i * ParamsetGamma[,'gamma4'])
    naturalspace <- cbind(Mod(temp[,1])^2 + 1, Arg(temp[,1])/2)
    if(whichLogged[2]-whichLogged[1]>1){
      Paramset <- cbind(exp(ParamsetGamma[, whichLogged[1]]), ParamsetGamma[,-c(whichLogged, whichAniso)], exp(ParamsetGamma[, whichLogged[2]]), naturalspace)
    }else{
      Paramset <- cbind(exp(ParamsetGamma[, whichLogged]), ParamsetGamma[,-c(whichLogged, whichAniso)], naturalspace)
    }
  } 
  
  colnames(Paramset) <- names(Mle)  
  toAdd = setdiff(c('range','shape','nugget','anisoRatio', 'anisoAngleRadians'), names(Mle))
  otherParams = matrix(Model$opt$mle[toAdd], nrow=nrow(Paramset), ncol = length(toAdd),
                       dimnames = list(rownames(Paramset), toAdd), byrow=TRUE)
  Params2 <- cbind(Paramset, otherParams)
  
  result2<-gpuRandom::getProfLogL(data= Model$data,
                                  formula=Model$model$formula,
                                  coordinates=Model$data@coords,
                                  params=Params2,
                                  boxcox = Model$parameters['boxcox'],
                                  type = "double",
                                  NparamPerIter=100,
                                  gpuElementsOnly=FALSE,
                                  reml=FALSE,
                                  Nglobal=c(128,64),
                                  Nlocal=c(16,16),
                                  NlocalCache=2000)
  
  
  A <- result2$LogLik[, 2]
  FirstDeri <- rep(0, length(Mle))
  names(FirstDeri) <- names(MleGamma)
  
  FirstDeri[1] <- (A[1] - A[2])/(2*deltas[1])
  FirstDeri[2] <- (A[3] - A[4])/(2*deltas[2])
  FirstDeri[3] <- (A[5] - A[6])/(2*deltas[3])
  FirstDeri[4] <- (A[7] - A[8])/(2*deltas[4])
  FirstDeri[5] <- (A[9] - A[10])/(2*deltas[5])
  
  
  Sigma <- -solve(HessianMat)
  index <- which(abs(FirstDeri) > 0.01)
  newMle <- Mle
  
  for(i in 1:length(index)){
    newMle[index[i]] <- Mle[index[i]] + Sigma[index[i],index[i]]*FirstDeri[index[i]]
  } 
  
  # for(i in 1:length(Mle)){
  #   newMle[i] <- Mle[i] + Sigma[i,i]*FirstDeri[i]
  # } 
  
  ## check if newMle's are on boundary
  if(newMle['nugget'] < 0){
    newMle['nugget'] = 0
  }
  
  if(newMle['shape'] < 0){
    newMle['shape'] = 0
  }
  
  # check mle shape
  if(newMle['shape'] >= kappa){
    newMle['shape'] = 1/newMle['shape']
    # if(Mle['shape'] <= delta)
    #   Mle['shape'] = delta
    #parToLog <- parToLog[!parToLog %in% 'shape']
  }
  
  
  output = list(Gradiant= FirstDeri,
                HessianMat = HessianMat,
                Sigma = Sigma,
                originalPoint = Mle,
                centralPoint = newMle,
                parToLog = parToLog,
                whichAniso = whichAniso,
                data = cbind(Params1, result1$LogLik[, 2]))
  
  output
}

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
#' @title set loglikelihood locations
#' @useDynLib gpuRandom
#' @export

    configParams <- function(Model,
                             alpha=c(0.0001, 0.01, 0.1, 0.2, 0.5, 0.8, 0.95, 0.99),
                             Mle = NULL# a vector of confidence levels 1-alpha
    ){
      

      
      ## get the Hessian
      output <- getHessian(Model = Model,
                           Mle = Mle)
   
      newMle <- output$centralPoint
      Sigma <- output$Sigma
      whichAniso <- output$whichAniso
      parToLog <- output$parToLog
      
      ## always find the stationary points/new Mles 
      # for(i in 1:length(index)){
      #   newMle[index[i]] <- Mle[index[i]] + Sigma[index[i],index[i]]*FirstDeri[index[i]]
      # } 
      if(newMle['nugget'] < 0.01){
        parToLog <- parToLog[!parToLog %in% 'nugget']
      }  
      whichLogged = which(names(newMle) %in% parToLog)
      
      
      if('anisoRatio' %in% names(newMle)){
        if(!('anisoAngleRadians' %in% names(newMle))){
          stop('anisoRatio and anisoAngleRadians must be together')
        }
        if(newMle['anisoRatio'] <= 1){
          gamma3 <-  unname(sqrt(1/newMle['anisoRatio']-1) * cos(2*newMle['anisoAngleRadians']))
          gamma4 <-  unname(sqrt(1/newMle['anisoRatio']-1) * sin(2*newMle['anisoAngleRadians']))
        }else{
          gamma3 <-  unname(sqrt(newMle['anisoRatio']-1) * cos(2*newMle['anisoAngleRadians']))
          gamma4 <-  unname(sqrt(newMle['anisoRatio']-1) * sin(2*newMle['anisoAngleRadians']))
        }
        aniso <- c(gamma3 = gamma3, gamma4 = gamma4)
        if(whichLogged[2]-whichLogged[1]>1){
          newMleGamma = c(log(newMle[whichLogged[1]]), newMle[-c(whichLogged, whichAniso)],log(newMle[whichLogged[2]]), aniso) 
        }else{
          newMleGamma = c(log(newMle[whichLogged]), newMle[-c(whichLogged, whichAniso)], aniso)
        }
      }else{
        newMleGamma = c(log(newMle[whichLogged]),  newMle[-whichLogged])
      }
      names(newMleGamma)[whichLogged] = paste("log(", names(newMle)[whichLogged], ")",sep="")
      
      

       
      ## fix first derivative ends  
      eig <- eigen(Sigma)
      out_list <- list()
      
      if(length(newMle)==2){
        pointsSphere = exp(1i*seq(0, 2*pi, len=25))
        pointsSphere = pointsSphere[-length(pointsSphere)]
        pointsSphere2d = cbind(Re(pointsSphere), Im(pointsSphere))
        #plot(pointsSphere2d)
        if('anisoRatio' %in% names(newMle)){
          for(i in 1:length(alpha)){
            clevel <- stats::qchisq(1 - alpha[i], df = 2)
            pointsEllipseGammaspace = t(sqrt(clevel) * eig$vectors %*% diag(sqrt(eig$values)) %*%  t(pointsSphere2d) + newMleGamma)
            colnames(pointsEllipseGammaspace) <- names(newMleGamma)
            temp <- as.data.frame(pointsEllipseGammaspace[,'gamma3'] + 1i * pointsEllipseGammaspace[,'gamma4'])
            naturalspace <- cbind(Mod(temp[,1])^2 + 1, Arg(temp[,1])/2)
            pointsEllipse <- cbind(exp(pointsEllipseGammaspace[,whichLogged]),naturalspace)
            colnames(pointsEllipse) <- names(newMle)
            out_list[[i]] = pointsEllipse
          }         
        }else{
          for(i in 1:length(alpha)){
            clevel <- stats::qchisq(1 - alpha[i], df = 2)
            pointsEllipseGammaspace = t(sqrt(clevel) * eig$vectors %*% diag(sqrt(eig$values)) %*%  t(pointsSphere2d) + newMleGamma)
            colnames(pointsEllipseGammaspace) <- names(newMleGamma)
            pointsEllipse <- cbind(exp(pointsEllipseGammaspace[,paste("log(", parToLog[whichLogged], ")",sep="")]))
            colnames(pointsEllipse) <- names(newMle)
            out_list[[i]] = pointsEllipse
          } 
        }
      }else if(length(newMle)==3){
        load('/home/ruoyong/diseasemapping/pkg/gpuRandom/data/coords3d.RData')
        if('anisoRatio' %in% names(newMle)){
          for(i in 1:length(alpha)){
            clevel <- stats::qchisq(1 - alpha[i], df = 3)
            pointsEllipseGammaspace = t(sqrt(clevel) * eig$vectors %*% diag(sqrt(eig$values)) %*%  t(coords3d) + newMleGamma)
            colnames(pointsEllipseGammaspace) <- names(newMleGamma)
            temp <- as.data.frame(pointsEllipseGammaspace[,'gamma3'] + 1i * pointsEllipseGammaspace[,'gamma4'])
            naturalspace <- cbind(Mod(temp[,1])^2 + 1, Arg(temp[,1])/2)
            pointsEllipse <- cbind(exp(pointsEllipseGammaspace[,whichLogged]),naturalspace)
            colnames(pointsEllipse) <- names(newMle)
            out_list[[i]] = pointsEllipse
          }   
        }else{
          for(i in 1:length(alpha)){
            clevel <- stats::qchisq(1 - alpha[i], df = 3)
            pointsEllipseGammaspace = t(sqrt(clevel) * eig$vectors %*% diag(sqrt(eig$values)) %*%  t(coords3d) + newMleGamma)
            colnames(pointsEllipseGammaspace) <- names(newMleGamma)
            pointsEllipse <- cbind(exp(pointsEllipseGammaspace[,paste("log(", parToLog[whichLogged], ")",sep="")]))
            colnames(pointsEllipse) <- names(newMle)
            out_list[[i]] = pointsEllipse
          }
        }
      }else if(length(newMle)==4){
        load('/home/ruoyong/diseasemapping/pkg/gpuRandom/data/coords4d.RData')
        for(i in 1:length(alpha)){
          clevel <- stats::qchisq(1 - alpha[i], df = 4)
          pointsEllipseGammaspace = t(sqrt(clevel) * eig$vectors %*% diag(sqrt(eig$values)) %*%  t(coords4d) + newMleGamma)
          colnames(pointsEllipseGammaspace) <- names(newMleGamma)
          temp <- pointsEllipseGammaspace[,'gamma3'] + 1i * pointsEllipseGammaspace[,'gamma4']
          pointsEllipse <- cbind(exp(pointsEllipseGammaspace[,whichLogged]), Mod(temp)^2 + 1, Arg(temp)/2)
          colnames(pointsEllipse) <- names(newMle)
          out_list[[i]] = pointsEllipse
        }
      }else if(length(newMle)==5){
        load('/home/ruoyong/diseasemapping/pkg/gpuRandom/data/coords5d.RData')
        # if(newMle['anisoRatio'] <= 1) {
        #    for(i in 1:length(alpha)){
        #       clevel <- stats::qchisq(1 - alpha[i], df = 5)
        #       pointsEllipseGammaspace = t(sqrt(clevel) * eig$vectors %*% diag(sqrt(eig$values)) %*%  t(coords5d) + newMleGamma)
        #       colnames(pointsEllipseGammaspace) <- names(newMleGamma)
        #       temp <- pointsEllipseGammaspace[,'gamma3'] + 1i * pointsEllipseGammaspace[,'gamma4']
        #       pointsEllipse <- cbind(exp(pointsEllipseGammaspace[,paste("log(", parToLog, ")",sep="")]), 1/(Mod(temp) + 1), Arg(temp)/2)
        #       colnames(pointsEllipse) <- names(newMle)
        #       out_list[[i]] = pointsEllipse
        #    }
        # }else{
        for(i in 1:length(alpha)){
          clevel <- stats::qchisq(1 - alpha[i], df = 5)
          pointsEllipseGammaspace = t(sqrt(clevel) * eig$vectors %*% diag(sqrt(eig$values)) %*%  t(coords5d) + newMleGamma)
          colnames(pointsEllipseGammaspace) <- names(newMleGamma)
          temp <- pointsEllipseGammaspace[,'gamma3'] + 1i * pointsEllipseGammaspace[,'gamma4']
          pointsEllipse <- cbind(exp(pointsEllipseGammaspace[,paste("log(", parToLog, ")",sep="")]), Mod(temp)^2 + 1, Arg(temp)/2)
          colnames(pointsEllipse) <- names(newMle)
          out_list[[i]] = pointsEllipse
        }
      }
      
      if(newMle['nugget'] == 0){
      for(i in 1:length(alpha)){
         vector <- out_list[[i]][,'nugget']
         for(j in 1:length(vector)/3){
           if(vector[j]<0){
             vector[j] = 0
           }
         }
         for(j in (length(vector)/3+1):length(vector)){
           if(vector[j]<0){
             vector[j] = stats::runif(1, 0, delta)
           }
         }
      }
      }
      
      if(abs(newMle['shape']) < 0.1){
      for(i in 1:length(alpha)){
        vector <- 1/out_list[[i]][,'shape']
        vector2 <- rep(1000, length(vector))
        out_list[[i]][,'shape'] = pmin(vector, vector2)
      }
      }

      names(out_list) <- paste0("alpha", alpha, sep="")
      out_list
    
   }
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
       
    
    
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    