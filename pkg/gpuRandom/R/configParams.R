#' @title set loglikelihood locations
#' @useDynLib gpuRandom
#' @export


    getHessian <- function(Model,
                           alpha=c(1e-6, 0.0005, 0.01, 0.1, 0.2, 0.5, 0.8, 0.95, 0.99),  # a vector of confidence levels 1-alpha
                           Mle = NULL){
  
  ## get the Hessian
  Theorder <- c('range', 'shape', 'nugget', 'anisoRatio', 'anisoAngleRadians')
  parToLog <- c("range", 'shape', "nugget")
  if(is.null(Mle)){
    Mle <- Model$optim$options[,'opt']
  }
  Mle <- Mle[order(match(names(Mle), Theorder))]
  Mle <- Mle[!names(Mle) %in% c('boxcox')]
  
  # if (!(identical(names(Mle), c('range', 'nugget')) |
  #     identical(names(Mle), c('range', 'nugget', 'anisoRatio', 'anisoAngleRadians'))|
  #     identical(names(Mle), c('range', 'nugget', 'shape', 'anisoRatio', 'anisoAngleRadians'))
  #     )){
  #    stop('Either 2 or 4 or 5 paramters to estimate')
  # }
  
  #parToLog = intersect(names(Mle), parToLog)
  whichLogged = which(names(Mle) %in% parToLog)
  whichAniso = which(names(Mle) %in% c('anisoRatio', 'anisoAngleRadians'))
  
  if('anisoRatio' %in% names(Mle)){
    if(Mle['anisoRatio'] <= 1) {
      gamma3 <-  unname(sqrt(1/Mle['anisoRatio']-1) * cos(2*Mle['anisoAngleRadians']))
      gamma4 <-  unname(sqrt(1/Mle['anisoRatio']-1) * sin(2*Mle['anisoAngleRadians']))
    }else{
      gamma3 <-  unname(sqrt(Mle['anisoRatio']-1) * cos(2*Mle['anisoAngleRadians']))
      gamma4 <-  unname(sqrt(Mle['anisoRatio']-1) * sin(2*Mle['anisoAngleRadians']))
    }
    aniso <- c(gamma3 = gamma3, gamma4 = gamma4)
    MleGamma = c(log(Mle[whichLogged]), Mle[-c(whichLogged, whichAniso)], aniso)
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
    Paramset <- cbind(exp(ParamsetGamma[, whichLogged]), naturalspace)
  } 
  
  colnames(Paramset) <- names(Mle)  
  toAdd = setdiff(c('range','shape','nugget','anisoRatio', 'anisoAngleRadians'), names(Mle))
  otherParams = matrix(Model$opt$mle[toAdd], nrow=nrow(Paramset), ncol = length(toAdd),
                       dimnames = list(rownames(Paramset), toAdd), byrow=TRUE)
  Params <- cbind(Paramset, otherParams)
  
  
  
  result<-gpuRandom::getProfLogL(data= Model$data,
                                 formula=Model$model$formula,
                                 coordinates=Model$data@coords,
                                 params=Params,
                                 boxcox = Model$parameters['boxcox'],
                                 type = "double",
                                 NparamPerIter=400,
                                 gpuElementsOnly=FALSE,
                                 reml=FALSE,
                                 Nglobal=c(128,64),
                                 Nlocal=c(16,16),
                                 NlocalCache=2800)
  
  
  #result$paramsRenew
  A <- result$LogLik[-1, 2]
  Origin <- result$LogLik[1, 2]
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
  
  HessianMat
}

    
    
    
    
#' @title set loglikelihood locations
#' @useDynLib gpuRandom
#' @export
       configParams <- function(Model,
                                alpha=c(1e-6, 0.0005, 0.01, 0.1, 0.2, 0.5, 0.8, 0.95, 0.99),
                                Mle = NULL# a vector of confidence levels 1-alpha
    ){
      
      ## get the First derivative
      
      # order the params
      Theorder <- c('range', 'shape', 'nugget', 'anisoRatio', 'anisoAngleRadians')
      parToLog <- c("range", 'shape', "nugget")
      if(is.null(Mle)){
        Mle <- Model$optim$options[,'opt']
      }
      Mle <- Mle[order(match(names(Mle), Theorder))]
      Mle <- Mle[!names(Mle) %in% c('boxcox')]
      
      # check if shape and nugget's mle is close to 0
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
        MleGamma = c(log(Mle[whichLogged]), Mle[-c(whichLogged, whichAniso)], aniso)
      }else{
        MleGamma = c(log(Mle[whichLogged]),  Mle[-whichLogged])
      }
      names(MleGamma)[whichLogged] = paste("log(", names(Mle)[whichLogged], ")",sep="")
      
      
      # frst detivatives
      derivGridDf <- rbind(c(1,0,0,0,0),
                           c(-1, 0,0,0,0),
                           c(0, 1, 0,0,0),
                           c(0,-1, 0,0,0),
                           c(0, 0, 1,0,0),
                           c(0, 0,-1,0,0),
                           c(0, 0, 0,1,0),
                           c(0, 0, 0,-1,0),
                           c(0, 0, 0, 0,1),
                           c(0, 0, 0, 0,-1))

      deltas = rep(0.01, length(Mle))
      # names(deltas) = names(MleGamma)
      # deltas['gamma4'] = 0.02
      ParamsetGamma <- matrix(MleGamma, nrow=nrow(derivGridDf), ncol=length(Mle), byrow=TRUE, dimnames = list(NULL, names(MleGamma))) + 
        derivGridDf %*% diag(deltas)
      
      if(!('anisoRatio' %in% names(Mle))){
        Paramset <- cbind(exp(ParamsetGamma[,paste("log(", names(Mle)[whichLogged], ")",sep="")]), ParamsetGamma[,-whichLogged])   
      }else{
        temp <- as.data.frame(ParamsetGamma[,'gamma3'] + 1i * ParamsetGamma[,'gamma4'])
        naturalspace <- cbind(Mod(temp[,1])^2 + 1, Arg(temp[,1])/2)
        Paramset <- cbind(exp(ParamsetGamma[, whichLogged]), naturalspace)
      } 
      
      colnames(Paramset) <- names(Mle)  
      toAdd = setdiff(c('range','shape','nugget','anisoRatio', 'anisoAngleRadians'), names(Mle))
      otherParams = matrix(Model$opt$mle[toAdd], nrow=nrow(Paramset), ncol = length(toAdd),
                           dimnames = list(rownames(Paramset), toAdd), byrow=TRUE)
      Params <- cbind(Paramset, otherParams)
      
      
      result<-gpuRandom::getProfLogL(data= Model$data,
                                     formula=Model$model$formula,
                                     coordinates=Model$data@coords,
                                     params=Params,
                                     boxcox = Model$parameters['boxcox'],
                                     type = "double",
                                     NparamPerIter=400,
                                     gpuElementsOnly=FALSE,
                                     reml=FALSE,
                                     Nglobal=c(128,64),
                                     Nlocal=c(16,16),
                                     NlocalCache=2800)
      
      
      A <- result$LogLik[, 2]
      FirstDeri <- rep(0, length(Mle))
      names(FirstDeri) <- names(MleGamma)
      
      FirstDeri[1] <- (A[1] - A[2])/(2*deltas[1])
      FirstDeri[2] <- (A[3] - A[4])/(2*deltas[2])
      FirstDeri[3] <- (A[5] - A[6])/(2*deltas[3])
      FirstDeri[4] <- (A[7] - A[8])/(2*deltas[4])
      FirstDeri[5] <- (A[9] - A[10])/(2*deltas[5])
      ## get the First derivative ends      
      
      
      ## find the stationary points/new Mles if first derivative not close to 0
      if(any(abs(FirstDeri) >= 0.01)){
        index <- which(abs(FirstDeri) >= 0.01)
        ToFix <- FirstDeri[index]
        
       for(i in 1:length(index)){
          theta0 <- Mle[index[i]]
          Htheta0 <- getHessian(Model,
                                     alpha=alpha,  # a vector of confidence levels 1-alpha
                                     Mle = theta0)
          ToFix[i] <- -(FirstDeri[index[i]] - Htheta0 * theta0)/Htheta0
          
       } 
          
          
        # }else if(length(ToFix)==2){
        # }
      
      Mle[index] <- ToFix
      }
       
      ## fix first derivative ends  
      
      
      # ## calculates Hessian for renewed Mles
      # if('shape' %in% names(Mle)){
      #   Mle1 <- Mle['shape']
      #   Mle2 <- Mle[-which(names(Mle)=='shape')]
      #   othersLength <- length(Mle2)
      #   
      #   shapeHessian <- getHessian(Model,
      #                              alpha=alpha,  # a vector of confidence levels 1-alpha
      #                              Mle = Mle1)
      #   otherHessian <- getHessian(Model,
      #                              alpha=alpha,  # a vector of confidence levels 1-alpha
      #                              Mle = Mle2)
      #   
      #   
      #   temp1 <- rbind(otherHessian[1,], rep(0,othersLength), otherHessian[2:othersLength,])
      #   
      #   HessianMat <- cbind(temp1[,1], rep(0,othersLength+1), temp1[,2:othersLength])
      #   
      #   HessianMat[2,2] <- shapeHessian
      #   
      #   colnames(HessianMat)[1:2] <- c('log(range)','log(shape)')
      #   rownames(HessianMat) <- colnames(HessianMat)
      # }else{
        HessianMat <- getHessian(Model,
                                 alpha=alpha,  # a vector of confidence levels 1-alpha
                                 Mle = Mle2)
      # }
      
      
      
      Sigma <- -solve(HessianMat)
      eig <- eigen(Sigma)
      out_list <- list()
      
      if(length(Mle)==2){
        pointsSphere = exp(1i*seq(0, 2*pi, len=25))
        pointsSphere = pointsSphere[-length(pointsSphere)]
        pointsSphere2d = cbind(Re(pointsSphere), Im(pointsSphere))
        #plot(pointsSphere2d)
        if('anisoRatio' %in% names(Mle)){
          for(i in 1:length(alpha)){
            clevel <- stats::qchisq(1 - alpha[i], df = 2)
            pointsEllipseGammaspace = t(sqrt(clevel) * eig$vectors %*% diag(sqrt(eig$values)) %*%  t(pointsSphere2d) + MleGamma)
            colnames(pointsEllipseGammaspace) <- names(MleGamma)
            temp <- as.data.frame(pointsEllipseGammaspace[,'gamma3'] + 1i * pointsEllipseGammaspace[,'gamma4'])
            naturalspace <- cbind(Mod(temp[,1])^2 + 1, Arg(temp[,1])/2)
            pointsEllipse <- cbind(exp(pointsEllipseGammaspace[,whichLogged]),naturalspace)
            colnames(pointsEllipse) <- names(Mle)
            out_list[[i]] = pointsEllipse
          }         
        }else{
          for(i in 1:length(alpha)){
            clevel <- stats::qchisq(1 - alpha[i], df = 2)
            pointsEllipseGammaspace = t(sqrt(clevel) * eig$vectors %*% diag(sqrt(eig$values)) %*%  t(pointsSphere2d) + MleGamma)
            colnames(pointsEllipseGammaspace) <- names(MleGamma)
            pointsEllipse <- cbind(exp(pointsEllipseGammaspace[,paste("log(", parToLog[whichLogged], ")",sep="")]))
            colnames(pointsEllipse) <- names(Mle)
            out_list[[i]] = pointsEllipse
          } 
        }
      }else if(length(Mle)==3){
        load('/home/ruoyong/diseasemapping/pkg/gpuRandom/data/coords3d.RData')
        if('anisoRatio' %in% names(Mle)){
          for(i in 1:length(alpha)){
            clevel <- stats::qchisq(1 - alpha[i], df = 3)
            pointsEllipseGammaspace = t(sqrt(clevel) * eig$vectors %*% diag(sqrt(eig$values)) %*%  t(coords3d) + MleGamma)
            colnames(pointsEllipseGammaspace) <- names(MleGamma)
            temp <- as.data.frame(pointsEllipseGammaspace[,'gamma3'] + 1i * pointsEllipseGammaspace[,'gamma4'])
            naturalspace <- cbind(Mod(temp[,1])^2 + 1, Arg(temp[,1])/2)
            pointsEllipse <- cbind(exp(pointsEllipseGammaspace[,whichLogged]),naturalspace)
            colnames(pointsEllipse) <- names(Mle)
            out_list[[i]] = pointsEllipse
          }   
        }else{
          for(i in 1:length(alpha)){
            clevel <- stats::qchisq(1 - alpha[i], df = 3)
            pointsEllipseGammaspace = t(sqrt(clevel) * eig$vectors %*% diag(sqrt(eig$values)) %*%  t(coords3d) + MleGamma)
            colnames(pointsEllipseGammaspace) <- names(MleGamma)
            pointsEllipse <- cbind(exp(pointsEllipseGammaspace[,paste("log(", parToLog[whichLogged], ")",sep="")]))
            colnames(pointsEllipse) <- names(Mle)
            out_list[[i]] = pointsEllipse
          }
        }
      }else if(length(Mle)==4){
        load('/home/ruoyong/diseasemapping/pkg/gpuRandom/data/coords4d.RData')
        for(i in 1:length(alpha)){
          clevel <- stats::qchisq(1 - alpha[i], df = 4)
          pointsEllipseGammaspace = t(sqrt(clevel) * eig$vectors %*% diag(sqrt(eig$values)) %*%  t(coords4d) + MleGamma)
          colnames(pointsEllipseGammaspace) <- names(MleGamma)
          temp <- pointsEllipseGammaspace[,'gamma3'] + 1i * pointsEllipseGammaspace[,'gamma4']
          pointsEllipse <- cbind(exp(pointsEllipseGammaspace[,whichLogged]), Mod(temp)^2 + 1, Arg(temp)/2)
          colnames(pointsEllipse) <- names(Mle)
          out_list[[i]] = pointsEllipse
        }
      }else if(length(Mle)==5){
        load('/home/ruoyong/diseasemapping/pkg/gpuRandom/data/coords5d.RData')
        # if(Mle['anisoRatio'] <= 1) {
        #    for(i in 1:length(alpha)){
        #       clevel <- stats::qchisq(1 - alpha[i], df = 5)
        #       pointsEllipseGammaspace = t(sqrt(clevel) * eig$vectors %*% diag(sqrt(eig$values)) %*%  t(coords5d) + MleGamma)
        #       colnames(pointsEllipseGammaspace) <- names(MleGamma)
        #       temp <- pointsEllipseGammaspace[,'gamma3'] + 1i * pointsEllipseGammaspace[,'gamma4']
        #       pointsEllipse <- cbind(exp(pointsEllipseGammaspace[,paste("log(", parToLog, ")",sep="")]), 1/(Mod(temp) + 1), Arg(temp)/2)
        #       colnames(pointsEllipse) <- names(Mle)
        #       out_list[[i]] = pointsEllipse
        #    }
        # }else{
        for(i in 1:length(alpha)){
          clevel <- stats::qchisq(1 - alpha[i], df = 5)
          pointsEllipseGammaspace = t(sqrt(clevel) * eig$vectors %*% diag(sqrt(eig$values)) %*%  t(coords5d) + MleGamma)
          colnames(pointsEllipseGammaspace) <- names(MleGamma)
          temp <- pointsEllipseGammaspace[,'gamma3'] + 1i * pointsEllipseGammaspace[,'gamma4']
          pointsEllipse <- cbind(exp(pointsEllipseGammaspace[,paste("log(", parToLog, ")",sep="")]), Mod(temp)^2 + 1, Arg(temp)/2)
          colnames(pointsEllipse) <- names(Mle)
          out_list[[i]] = pointsEllipse
        }
        
      }
      
      #out_list[length(alpha)+1] = HessianMat
      names(out_list) <- paste0("alpha", alpha, sep="")
      out_list
    
   }
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
       
    
    
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    