library('geostatsp')
data("swissRain")
swissRain$elevation = extract(swissAltitude, swissRain)
swissRes =  lgm( formula=rain~ elevation, 
                     data=swissRain, 
                     grid=20,
                     covariates=swissAltitude, 
                     reml = FALSE,
                     fixBoxcox=FALSE, fixShape=FALSE, fixNugget = FALSE,  #Set to FALSE to estimate the nugget effect parameter.
                     aniso=TRUE )
swissRes$summary[,c('estimate','ci0.025', 'ci0.975')]
swissRes$summary[,c('estimate','ci0.1', 'ci0.9')]


## ---gpuRandom------range vs nugget contour--------------------------------------------------------------
## set params
library(gpuRandom)
# params = expand.grid(
#   range = c(exp(seq(log(15), log(swissRes$parameters['range']/1000), len=6))*1000, exp(seq(log(swissRes$parameters['range']/1000+5), log(250), len=6))*1000),
#   shape = c(seq(0.1, swissRes$parameters['shape'], len=5), seq(swissRes$parameters['shape']+0.5, 5, len=5)),
#   anisoRatio =c( seq(1, swissRes$parameters['anisoRatio'], len=5), seq(swissRes$parameters['anisoRatio']+1, 12, len=5) ),
#   anisoAngleDegrees = c( seq(27,swissRes$parameters['anisoAngleDegrees'], len=5),  seq(swissRes$parameters['anisoAngleDegrees']+2, 46, len=5)),
#   nugget = c(seq(0, swissRes$optim$mle['nugget'], len=3), seq(swissRes$optim$mle['nugget']+0.5, 8, len=8))
# )

ParamList = list(
  range = c(exp(seq(log(15), log(swissRes$parameters['range']/1000), len=6))*1000, exp(seq(log(swissRes$parameters['range']/1000+5), log(240), len=5))*1000),
  shape = c(seq(0.1, swissRes$parameters['shape'], len=5), seq(swissRes$parameters['shape']+0.5, 4.5, len=4)),
  anisoRatio =c( seq(1, swissRes$parameters['anisoRatio'], len=5), seq(swissRes$parameters['anisoRatio']+1, 16, len=6) ),
  anisoAngleDegrees = c( seq(27,swissRes$parameters['anisoAngleDegrees'], len=5),  seq(swissRes$parameters['anisoAngleDegrees']+2, 46, len=5)),
  nugget = c(seq(0, swissRes$optim$mle['nugget'], len=3), seq(swissRes$optim$mle['nugget']+0.5, 8, len=8))
) 
params = do.call(expand.grid, ParamList)






output <- gpuRandom::likfitLgmCov(spatialmodel,
                                  params=params, # CPU matrix for now, users need to provide proper parameters given their specific need
                                  boxcox=c(seq(0, swissResl$parameters['boxcox'], len=5), seq(0.6, 1.5, len=3)),  # boxcox is always estimated
                                  paramToEstimate = c("range", "shape", "nugget", "sdNugget","anisoRatio", "anisoAngleDegrees", "boxcox"), #variance and regression parameters are always estimated if not given,
                                  cilevel=0.9,  # decimal
                                  type = "double",
                                  reml=FALSE, 
                                  NparamPerIter=400,
                                  Nglobal=c(128,64),
                                  Nlocal=c(16,16),
                                  NlocalCache=2800,
                                  verbose=FALSE)
output$summary

result1<-gpuRandom::getProfLogL(data=swissRain,
                                formula=rain~ elevation,
                                coordinates=swissRain@coords,
                                params=params, 
                                boxcox = c(seq(0, swissRes$parameters['boxcox'], len=4), seq(0.6, 1.5, len=2)),
                                type = "double",
                                NparamPerIter=400,
                                gpuElementsOnly = FALSE,
                                reml=FALSE, 
                                Nglobal=c(128,64),
                                Nlocal=c(16,16),
                                NlocalCache=2800)



result1$Infindex
predictors = result1$predictors
LogLik = result1$LogLik # cpu matrix
XVYXVX = result1$XVYXVX  # cpu matrix
ssqResidual = result1$ssqResidual  # cpu matrix
paramToEstimate = c("range","nugget", "boxcox")
cilevel=0.9  # decimal
params = result1$params # cpu matrix,
boxcox = result1$boxcox  # boxcox vallues, consistent with other functions
Ndata = result1$Ndata
Nobs = result1$Nobs
Ncov = result1$Ncov
reml=FALSE
verbose=FALSE


maximum <- max(LogLik)
breaks = maximum - qchisq(cilevel,  df = 1)/2
index <- which(LogLik == max(LogLik, na.rm = TRUE), arr.ind = TRUE)
#################sigma hat#########################
if(reml==FALSE)  {
  Table["sdSpatial",1] <- sqrt(ssqResidual[index[1],index[2]]/Nobs)
}else{         
  Table["sdSpatial",1] <- sqrt(ssqResidual[index[1],index[2]]/(Nobs - Ncov))
}

params <- cbind(sqrt(params[,"nugget"]) * Table["sdSpatial",1], params)
colnames(params)[1] <- 'sdNugget'



par(mfrow = c(2, 2))
par(mar = c(2.5, 3.5, 2, 0.5))
par(mgp = c(1.5, 0.5, 0))
par(oma = c(3, 0, 3, 0))
result = data.table::as.data.table(cbind(LogLik, params[,"range"]))
colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 3) ,sep = ''), "range")
profileLogLik <- result[, .(profile=max(.SD)), by=range]
f1 <- approxfun(profileLogLik$range/1000, profileLogLik$profile-breaks)  
plot(profileLogLik$range/1000, profileLogLik$profile-breaks, ylab= "proLogL", xlab='range/1000', xlim=c(15, 240), xaxt='n')
curve(f1(x), add = TRUE, col = 2, n = 1001) 
axis(1,at=seq(15, 240, len=4),labels=TRUE)
abline(h =0, lty = 2)
lower = min(profileLogLik$range/1000);upper = max(profileLogLik$range/1000)
ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
abline(v =ci[1], lty = 2)
abline(v =ci[2], lty = 2)


result = as.data.table(cbind(LogLik, params[,"shape"]))
colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 3) ,sep = ''), "shape")
profileLogLik <- result[, .(profile=max(.SD)), by=shape]
plot(profileLogLik$shape, profileLogLik$profile-breaks,  ylab= "proLogL", xlab="shape")
f1 <- approxfun(profileLogLik$shape, profileLogLik$profile-breaks) 
curve(f1, add = TRUE, col = 2) 
abline(h =0, lty = 2)
lower = min(profileLogLik$shape);upper = max(profileLogLik$shape)
ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
abline(v =ci[1], lty = 2)



result = as.data.table(cbind(LogLik, params[,"sdNugget"]))
colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 2) ,sep = ''), "sdNugget")
profileLogLik <-result[, .(profile=max(.SD)), by=sdNugget]
plot(profileLogLik$sdNugget, profileLogLik$profile-breaks,  ylab= "proLogL", xlab=bquote(tau), xlim=c(0, 8), xaxt='n')
f1 <- approxfun(profileLogLik$sdNugget, profileLogLik$profile-breaks)  
curve(f1, add = TRUE, col = 2, n=1001) 
axis(1,at=seq(0, 8, len=5),labels=TRUE)
abline(h =0, lty = 2)
lower = min(profileLogLik$sdNugget)
upper = max(profileLogLik$sdNugget)
ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
abline(v =ci[1], lty = 2)




result = as.data.table(cbind(LogLik, params[,"nugget"]))
colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 2) ,sep = ''), "nugget")
profileLogLik <-result[, .(profile=max(.SD)), by=nugget]
plot(profileLogLik$nugget, profileLogLik$profile-breaks, ylab= "proLogL", xlab=expression(nu^2))
f1 <- approxfun(profileLogLik$nugget, profileLogLik$profile-breaks)  
curve(f1, add = TRUE, col = 2, n=1001) 
abline(h =0, lty = 2)
lower = min(profileLogLik$nugget)
upper = max(profileLogLik$nugget)
MLE <-  optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum
ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
abline(v =ci[1], lty = 2)







  result = as.data.table(cbind(LogLik, params[,"anisoRatio"]))
  colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 2) ,sep = ''), "anisoRatio")
  profileLogLik <- result[, .(profile=max(.SD)), by=anisoRatio]
  plot(profileLogLik$anisoRatio, profileLogLik$profile-breaks,ylab= "proLogL", xlab='anisoRatio')
  f1 <- approxfun(profileLogLik$anisoRatio, profileLogLik$profile-breaks)  
  curve(f1, add = TRUE, col = 2, n=1001) 
  abline(h =0, lty = 2)
  lower = min(profileLogLik$anisoRatio)
  upper = max(profileLogLik$anisoRatio)
  ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
  abline(v =ci[1], lty = 2)
  abline(v =ci[2], lty = 2)

  
  
  
  
  
if('anisoAngleRadians' %in% paramToEstimate){
  result = as.data.table(cbind(LogLik, params[,"anisoAngleRadians"]))
  colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 2) ,sep = ''), "anisoAngleRadians")
  profileLogLik <- result[, .(profile=max(.SD)), by=anisoAngleRadians]
  plot(profileLogLik$anisoAngleRadians, profileLogLik$profile-breaks, main="Profile LogL, y axis adjusted", ylab= "proLogL-breaks", xlab='anisoAngleRadians')
  f1 <- approxfun(profileLogLik$anisoAngleRadians, profileLogLik$profile-breaks)
  # f1 <- splinefun(profileLogLik$anisoAngleRadians, profileLogLik$profile-breaks, method = "monoH.FC")
  curve(f1, add = TRUE, col = 2, n=1001) 
  abline(h =0, lty = 2)
  #myplots[['anisoAngleRadians']] <- plot.Degrees
  # anisoAngleRadiansresults <- optim(0.1, f1, method = "L-BFGS-B",lower = 0.1, upper = 1.5, hessian = FALSE, control=list(fnscale=-1) )
  lower = min(profileLogLik$anisoAngleRadians)
  upper = max(profileLogLik$anisoAngleRadians)
  MLE <- round(optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum,digits=4)
  # maxvalue <- anisoAngleRadiansresults$objective
  # breaks = maxvalue - qchisq(cilevel,  df = 1)/2
  ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
  if(length(ci)==1){
    if( ci > MLE){
      ci <- c(lower, ci)
      message("did not find lower ci")
    }else{
      ci <- c(ci, upper)
      message("did not find upper ci")}
  }
  
  if(length(ci)==0 | length(ci)>2){
    warning("require a better param matrix")
    ci <- c(NA, NA)
  }
  Table["anisoAngleRadians",] <- c(MLE, ci)
  
  
}







if('anisoAngleDegrees' %in% paramToEstimate){
  result = as.data.table(cbind(LogLik, params[,"anisoAngleDegrees"]))
  colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 2) ,sep = ''), "anisoAngleDegrees")
  profileLogLik <- result[, .(profile=max(.SD)), by=anisoAngleDegrees]
  plot(profileLogLik$anisoAngleDegrees, profileLogLik$profile-breaks, main="Profile LogL, y axis adjusted", ylab= "proLogL-breaks", xlab='anisoAngleDegrees')
  f1 <- approxfun(profileLogLik$anisoAngleDegrees, profileLogLik$profile-breaks)
  # f1 <- splinefun(profileLogLik$anisoAngleDegrees, profileLogLik$profile-breaks, method = "monoH.FC")
  curve(f1, add = TRUE, col = 2, n=1001) 
  abline(h =0, lty = 2)
  #myplots[['anisoAngleDegrees']] <- plot.Degrees
  # anisoAngleDegreesresults <- optim(0.1, f1, method = "L-BFGS-B",lower = 0.1, upper = 1.5, hessian = FALSE, control=list(fnscale=-1) )
  lower = min(profileLogLik$anisoAngleDegrees)
  upper = max(profileLogLik$anisoAngleDegrees)
  MLE <- round(optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum,digits=4)
  # maxvalue <- anisoAngleDegreesresults$objective
  # breaks = maxvalue - qchisq(cilevel,  df = 1)/2
  ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
  if(length(ci)==1){
    if( ci > MLE){
      ci <- c(lower, ci)
      message("did not find lower ci")
    }else{
      ci <- c(ci, upper)
      message("did not find upper ci")}
  }
  
  if(length(ci)==0 | length(ci)>2){
    warning("require a better param matrix")
    ci <- c(NA, NA)
  }
  Table["anisoAngleDegrees",] <- c(MLE, ci)
  
  
}



result = as.data.table(cbind(LogLik, params[,"anisoRatio"]))
colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 2) ,sep = ''), "anisoRatio")
profileLogLik <- result[, .(profile=max(.SD)), by=anisoRatio]
plot(profileLogLik$anisoRatio, profileLogLik$profile-breaks, ylab= "proLogL", xlab='anisoRatio')
f1 <- approxfun(profileLogLik$anisoRatio, profileLogLik$profile-breaks)  
curve(f1, add = TRUE, col = 2, n=1001) 
abline(h =0, lty = 2)
lower = min(profileLogLik$anisoRatio)
upper = max(profileLogLik$anisoRatio)
MLE <- round(optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum,digits=4)
ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)























