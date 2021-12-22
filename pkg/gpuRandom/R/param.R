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


  dynamicParams <- function(Summary,      #lgm model
                            optimMle,
                            paramToEstimate){
  
  
    # temp = rbind(c('range', Summary['range','Estimated']),
    #       c('shape', Summary['shape','Estimated']),
    #       c('nugget', Summary['sdNugget','Estimated']),
    #       c('anisoRatio', Summary['anisoRatio','Estimated']),
    #       c('anisoAngleRadians', Summary['anisoAngleRadians','Estimated']),
    #       c('anisoAngleDegrees', Summary['anisoAngleDegrees','Estimated']))
    # 
    # paramToEstimate = temp[which(temp[,2]=='TRUE')]
    
    
    ParamList = list(
      range = c(seq(Summary["range", 'ci0.005'], Summary["range", 'estimate'], len = 10),   seq(Summary["range", 'estimate']*1.1, Summary["range", 'ci0.995'], len=10)),
      nugget = c(seq(0, optimMle['nugget'], len=5),  seq(optimMle['nugget']*1.2, Summary["sdNugget", 'ci0.995'], len=5)),
      shape = c(seq(Summary["shape", 'ci0.005'], Summary["shape", 'estimate'], len = 10),   seq(Summary["shape", 'estimate']*1.1, Summary["shape", 'ci0.995'], len=10)),
      anisoRatio = c(seq(Summary['anisoRatio', 'ci0.005'], Summary['anisoRatio', 'estimate'], len = 10),   seq(Summary['anisoRatio', 'estimate']*1.1, Summary['anisoRatio', 'ci0.995'], len=10)),
      anisoAngleRadians = c(seq(Summary['anisoAngleRadians', 'ci0.005'], Summary['anisoAngleRadians', 'estimate'], len = 10),   seq(Summary['anisoAngleRadians', 'estimate']*1.1, Summary['anisoAngleRadians', 'ci0.995'], len=10)),
      anisoAngleDegrees = c(seq(Summary['anisoAngleDegrees', 'ci0.005'], Summary['anisoAngleDegrees', 'estimate'], len = 10),   seq(Summary['anisoAngleDegrees', 'estimate']*1.1, Summary['anisoAngleDegrees', 'ci0.995'], len=10))
    )
    
    ParamList = ParamList[paramToEstimate]
    
    params = do.call(expand.grid, ParamList)
    
    boxcox = c(seq(Summary["boxcox", 'ci0.005'], Summary["boxcox", 'estimate'], len=3), seq(Summary["boxcox", 'estimate']+0.3, Summary["boxcox", 'ci0.995'], len=4))
    

  
    output = list(params = params,
                  boxcox = boxcox,
                  ParamList = ParamList)

    output
    
  #theOrder = c('range','shape','variance','nugget','anisoRatio','anisoAngleRadians')
  
}













