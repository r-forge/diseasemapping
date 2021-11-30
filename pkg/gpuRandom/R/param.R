library('geostatsp')
data("swissRain")
swissRain$sqrtY = sqrt(swissRain$rain)
swissRain$elevation = extract(swissAltitude, swissRain)

swissRes =  lgm( formula=sqrtY~ elevation, 
                 data=swissRain, 
                 grid=20,#covariates=swissAltitude, 
                 reml = FALSE,
                 fixBoxcox=FALSE, fixShape=FALSE, fixNugget = FALSE,  #Set to FALSE to estimate the nugget effect parameter.
                 aniso=TRUE )


swissRes$summary[,c('estimate','ci0.005', 'ci0.995')]

modSummary = swissRes$summary

dynamicParams <- function(modSummary
                          ){
  
  params = expand.grid(
           range = seq(modSummary["range", 'ci0.005'], modSummary["range", 'ci0.995'], len = 20),
           shape = seq(modSummary["shape", 'ci0.005'], modSummary["shape", 'ci0.995'], len = 20),
           variance = 1,
           nugget = seq(modSummary["sdNugget", 'ci0.005'], modSummary["sdNugget", 'ci0.995'], len = 10),
           anisoRatio = seq(modSummary["anisoRatio", 'ci0.005'], modSummary["anisoRatio", 'ci0.995'], len = 15),
           anisoAngleRadians = seq(modSummary["anisoAngleRadians", 'ci0.005'], modSummary["anisoAngleRadians", 'ci0.995'], len = 15)
  )
  
  
  
  
  
  #theOrder = c('range','shape','variance','nugget','anisoRatio','anisoAngleRadians')
  
  
  
  
  
  

  
  
  
}













