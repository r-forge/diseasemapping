#' @title Estimate profile Log-likelihood for covariance parameters and lambda
#' @useDynLib gpuRandom
#' @export



# library('geostatsp')
# data("swissRain")
# swissRain$sqrtY = sqrt(swissRain$rain)
# swissRain$elevation = extract(swissAltitude, swissRain)
# 
# swissRes =  lgm( formula=sqrtY~ elevation,
#                  data=swissRain,
#                  grid=20,#covariates=swissAltitude,
#                  reml = FALSE,
#                  fixBoxcox=FALSE, fixShape=FALSE, fixNugget = FALSE,  #Set to FALSE to estimate the nugget effect parameter.
#                  aniso=TRUE )
# 
# 
# swissRes$summary[,c('estimate','ci0.005', 'ci0.995')]
# Summary = swissRes$summary


  ParamsFromLgm <- function(Summary,      #lgm model
                            optimMle,
                            paramToEstimate){
  
  

    
    if(is.element('anisoAngleRadians',paramToEstimate)){
      paramToEstimate = union(paramToEstimate, 'anisoAngleDegrees')
    }
    
    if(is.element('anisoAngleDegrees',paramToEstimate)){
      paramToEstimate = union(paramToEstimate, 'anisoAngleRadians')
    }
    
    
    
    optimMle_degree <-  Summary['anisoAngleDegrees', 'estimate']
    
    ParamList1 = list(
      range = c(seq(Summary["range", 'ci0.005'], optimMle['range'], len = 10),   seq(optimMle['range']*1.1, Summary["range", 'ci0.995'], len=10)),
      shape = c(seq(0.1, optimMle['shape'], len = 8),   seq(optimMle['shape']*1.1, Summary["shape", 'ci0.995'], len=8)),
      nugget = c(seq(0, optimMle['nugget'], len=5),    seq(optimMle['nugget']*1.2, Summary["sdNugget", 'ci0.995'], len=5)),
      anisoRatio = c(seq(Summary['anisoRatio', 'ci0.005'], optimMle['anisoRatio'], len = 10),   seq(optimMle['anisoRatio']*1.1, Summary['anisoRatio', 'ci0.995'], len=10)),
      anisoAngleRadians = c(seq(Summary['anisoAngleRadians', 'ci0.005'], optimMle['anisoAngleRadians'], len = 10),   seq(optimMle['anisoAngleRadians']*1.1, Summary['anisoAngleRadians', 'ci0.995'], len=10)),
      anisoAngleDegrees = c(seq(Summary['anisoAngleDegrees', 'ci0.005'], optimMle_degree, len = 10),   seq(optimMle_degree*1.1, Summary['anisoAngleDegrees', 'ci0.995'], len=10))
    )
    
    ParamList2 = list(
      range = unname(optimMle['range']),
      shape = unname(optimMle['shape']), 
      nugget = unname(optimMle['nugget']),
      anisoRatio = unname(optimMle['anisoRatio']), 
      anisoAngleRadians = unname(optimMle['anisoAngleRadians']), 
      anisoAngleDegrees = optimMle_degree
    )
    
    theOrder = c('range','shape','nugget','anisoRatio','anisoAngleRadians', 'anisoAngleDegrees')
    noEstimate = setdiff(theOrder, paramToEstimate)
    
    ParamList = c(ParamList1[paramToEstimate], ParamList2[noEstimate])
    
    params = do.call(expand.grid, ParamList)
    
    boxcox = c(seq(Summary["boxcox", 'ci0.005'], Summary["boxcox", 'estimate'], len=3), seq(Summary["boxcox", 'estimate']+0.3, Summary["boxcox", 'ci0.995'], len=4))
    

  
    output = list(params = params,
                  boxcox = boxcox,
                  ParamList = ParamList)

    output
    

  
}



  
  #theOrder = c('range','shape','variance','nugget','anisoRatio','anisoAngleRadians')
  
  # temp = rbind(c('range', Summary['range','Estimated']),
  #       c('shape', Summary['shape','Estimated']),
  #       c('nugget', Summary['sdNugget','Estimated']),
  #       c('anisoRatio', Summary['anisoRatio','Estimated']),
  #       c('anisoAngleRadians', Summary['anisoAngleRadians','Estimated']),
  #       c('anisoAngleDegrees', Summary['anisoAngleDegrees','Estimated']))
  # 
  # paramToEstimate = temp[which(temp[,2]=='TRUE')]









