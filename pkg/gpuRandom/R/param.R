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
  
    if(isTRUE(paramToEstimate==NULL)){
      temp = rbind(c('range', Summary['range','Estimated']),
                   c('shape', Summary['shape','Estimated']),
                   c('nugget', Summary['sdNugget','Estimated']),
                   c('anisoRatio', Summary['anisoRatio','Estimated']),
                   c('anisoAngleRadians', Summary['anisoAngleRadians','Estimated']))
                   #c('anisoAngleDegrees', Summary['anisoAngleDegrees','Estimated']))
      
      paramToEstimate = temp[which(temp[,2]=='TRUE')]
    }
    
    paramToEstimate <- paramToEstimate[! paramToEstimate %in% 'boxcox']
    optimMle_degree <-  Summary['anisoAngleDegrees', 'estimate']
    
    ParamList1 = list(
      range = c(seq(Summary["range", 'ci0.005'], Summary["range", 'estimate'], len = 7),   seq(Summary["range", 'estimate']*1.1, Summary["range", 'ci0.995']*1.1, len=7))*1000,
      shape = c(seq(0.1, optimMle['shape'], len = 5),   seq(optimMle['shape']*1.1, Summary["shape", 'ci0.995'], len=5)),
      nugget = c(seq(0, optimMle['nugget'], len=5),    seq(optimMle['nugget']*1.2, optimMle['nugget']*4, len=5)),
      anisoRatio = c(seq(Summary['anisoRatio', 'ci0.005']/2, optimMle['anisoRatio'], len = 5),   seq(optimMle['anisoRatio']*1.1, Summary['anisoRatio', 'ci0.995'], len=5)),
      anisoAngleRadians = c(seq(Summary['anisoAngleRadians', 'ci0.005']/2, optimMle['anisoAngleRadians'], len = 6),   seq(optimMle['anisoAngleRadians']*1.1, Summary['anisoAngleRadians', 'ci0.995'], len=6)),
      anisoAngleDegrees = c(seq(Summary['anisoAngleDegrees', 'ci0.005'], optimMle_degree, len = 5),   seq(optimMle_degree+4, Summary['anisoAngleDegrees', 'ci0.995'], len=5))
    )
    
    ParamList2 = list(
      range = unname(optimMle['range']),
      shape = unname(optimMle['shape']), 
      nugget = unname(optimMle['nugget']),
      anisoRatio = unname(optimMle['anisoRatio']), 
      anisoAngleRadians = unname(optimMle['anisoAngleRadians']),
      anisoAngleDegrees = optimMle_degree
    )

    if(('anisoAngleRadians' %in% paramToEstimate) | ('anisoAngleDegrees' %in% paramToEstimate)){
      paramToEstimate2 = union(paramToEstimate, c('anisoAngleRadians','anisoAngleDegrees'))
    }

    
    theOrder = c('range','shape','nugget','anisoRatio','anisoAngleRadians','anisoAngleDegrees')
    noEstimate = setdiff(theOrder, paramToEstimate2)
    
    if('anisoAngleRadians' %in% paramToEstimate)
      paramToEstimate = paramToEstimate[! paramToEstimate %in% 'anisoAngleDegrees'] 

    
    if('anisoAngleDegrees' %in% paramToEstimate)
      paramToEstimate = paramToEstimate[! paramToEstimate %in% 'anisoAngleRadians'] 

    
    
    ParamList = c(ParamList1[paramToEstimate], ParamList2[noEstimate])
    
    params = do.call(expand.grid, ParamList)
    
    boxcox = c(seq(Summary["boxcox", 'ci0.005'], Summary["boxcox", 'estimate'], len=3), seq(Summary["boxcox", 'estimate']+0.3, Summary["boxcox", 'ci0.995'], len=3))
    

  
    output = list(params = params,
                  boxcox = boxcox,
                  ParamList = ParamList)

    output
    

  
}



  
  #theOrder = c('range','shape','variance','nugget','anisoRatio','anisoAngleRadians')
  










