maternCholSolve = function(param, obsCov, coordinates){
  

  if(TRUE){ # do this in R
#  covMat = geostatsp::matern(x=coordinates, param=param)	
#  cholCovMat = Matrix::chol(covMat)
    
    cholCovMat = geostatsp::matern(x=coordinates, 
        param=param, type='cholesky')
    
  # cholCovMat %*% t(cholCovMat) = covMat
  
  # cholCovInvX = cholCovMat^{-1} %*% covariates
  cholCovInvXY = Matrix::solve(cholCovMat, obsCov)
#  cholCovInvX = Matrix::solve(cholCovMat, covariates)
  cholCovInvX = cholCovInvXY[,-1]
  # cholCovInvY = cholCovMat^{-1} %*% observations
  #cholCovInvY = Matrix::solve(cholCovMat, observations)
  cholCovInvY = cholCovInvXY[,1]
  
  cholCovInvXcross = Matrix::crossprod(cholCovInvX)
  cholCovInvXcrossInv =       Matrix::solve(cholCovInvXcross)
  detCholCovInvXcross = Matrix::determinant(cholCovInvXcross)$modulus
  
  betaHat = as.vector(
      cholCovInvXcrossInv %*% 
          Matrix::crossprod(cholCovInvX, cholCovInvY)
  ) 
  
  resids = obsCov[,1] - as.vector(obsCov[,-1] %*% betaHat)
  # sigsqhat = resids' %*% Vinv %*% residsx
  #    =   resids' Linv' Linv resids
  cholCovInvResid = Matrix::solve(cholCovMat, resids)
  detCholCovMat = determinant(cholCovMat)$modulus
  
  Nobs = nrow(obsCov)
  Ncov = ncol(obsCov)-1
  Nadj = c(ml=Nobs, reml=Nobs-Ncov)
  totalSsq = as.vector(Matrix::crossprod(cholCovInvResid))
  totalVarHat = totalSsq/Nadj
  
  minusTwoLogLik =  
      Nadj * log(2*pi) +
      2*detCholCovMat  

  minusTwoLogLik = cbind(
      fixedVariance = minusTwoLogLik+totalSsq,
      estimatedVariance = 
           minusTwoLogLik + Nadj + Nadj * log(totalVarHat)
  )
  
  minusTwoLogLik['reml',] = minusTwoLogLik['reml',] + 2*detCholCovInvXcross
  
  varMle = outer(totalVarHat, param[c("variance","nugget")])
  
  names(betaHat) = colnames(cholCovInvXY)[-1]
  cholCovInvXcrossInv=as.matrix(cholCovInvXcrossInv)
  dimnames(cholCovInvXcrossInv)=list(names(betaHat),names(betaHat))
  
  varBetaHat = array(NA, 
      c(dim(cholCovInvXcrossInv),2,2),
      dimnames = c(dimnames(cholCovInvXcrossInv),
          list(
              c('fixedVariance','estimatedVariance'), 
              c('ml','reml')
          )
      )
  )
  
  varBetaHat[,,'fixedVariance','ml'] =
      varBetaHat[,,'fixedVariance','reml'] =
      cholCovInvXcrossInv
  
  varBetaHat[,,'estimatedVariance','ml'] =
      cholCovInvXcrossInv*totalVarHat['ml']
  varBetaHat[,,'estimatedVariance','reml'] =
      cholCovInvXcrossInv*totalVarHat['reml']
      
  resultR = list(
          minusTwoLogLik = minusTwoLogLik,
          varMle = varMle,
          betaHat = betaHat,
          varBetaHat = varBetaHat,
          totalVarHat = totalVarHat,
          resids=resids,
          detCholCovMat=detCholCovMat
          )
  } else { # don't use R
    resultR = NULL
  }
  if(FALSE){ # use C
    Nobs = nrow(obCov)
    Ncov = ncol(obsCov)-1
    resultC = #.C("maternCholSolve",
        list(
        obsCov = as.double(obsCov),
        Nobs = as.integer(Nobs),
        Ncov = as.integer(Ncov),
        coordinates = as.double(coordinates),
        Ncoordinates = as.integer(length(coordinates)),
        coordsDist = as.integer(class(coordinates)=='dist'),
        nugget=as.double(param['nugget']),
        haveNugget = as.integer(param['nugget']>0),
        range=as.double(param["range"]),
        shape=as.double(param["shape"]),
        variance=as.double(param["variance"]),
        anisoRatio = as.double(param["anisoRatio"]),
        anisoAngle = as.double(param["anisoAngleRadians"]),
        betaHat = as.double(rep(0, Ncov)),
        cholCovInvXcrossInv = as.double(rep(0, Ncov*Ncov)),
        totalSsq = as.double(0),
        detCholCovMat= as.double(0),
        detCholCovInvXcross= as.double(0)
      )
    resultC$cholCovInvXcrossInv = matrix(
          resultC$cholCovInvXcrossInv, Ncov, Ncov
          )
    resultC$resids = resultC$obsCov[1:Nobs]
    
    resultC$diff = NULL
    for(D in names(resultR)){
      resultC$diff = cbind(
          resultC$diff,
          range(resultC[[D]] - resultR[[D]])
      )
    }
    colnames(resultC$diff) = names(resultR)
  } else {# don't use C
    resultC = NULL
  }
  
  return(list(R=resultR, C=resultC))
}