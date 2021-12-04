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
# modSummary = swissRes$summary


  dynamicParams <- function(modSummary
                          ){
  
  params = expand.grid(
           range = seq(modSummary["range", 'ci0.005'], modSummary["range", 'ci0.995'], len = 20),
           shape = seq(modSummary["shape", 'ci0.025'], modSummary["shape", 'ci0.975'], len = 15),
           nugget = seq(modSummary["sdNugget", 'ci0.005'], modSummary["sdNugget", 'ci0.995'], len = 10),
           anisoRatio = seq(modSummary["anisoRatio", 'ci0.005']/5, modSummary["anisoRatio", 'ci0.995'], len = 10),
           anisoAngleDegrees = seq(modSummary["anisoAngleDegrees", 'ci0.005'], modSummary["anisoAngleDegrees", 'ci0.995'], len = 10)
  )
  
  
  params
  
  
  #theOrder = c('range','shape','variance','nugget','anisoRatio','anisoAngleRadians')
  
}













