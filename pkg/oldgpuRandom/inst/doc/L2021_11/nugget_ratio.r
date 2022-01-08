########################################################################################################################################
## ---gpuRandom------nugget vs ratio contour--------------------------------------------------------------
## set params
newParamListsn = list(
  nugget = c(seq(0,4,len=8), seq(5,20, len=8)), 
  anisoRatio = c(seq(0.1, 1, len=2), seq(2,9, len=8)),
  anisoAngleRadians = 0,
  shape=2,
  nugget = 4,  
  range = 1000
) 
paramssn = do.call(expand.grid, newParamListsn)

count = 0
breaksf = rep(0, 10)
makedata <- rep(0, 4)
estimates <- rep(0, 6)
## calculate logl
for (i in 1:15){
  result<-gpuRandom::likfitLgmCov2d(
    data= mydat,
    formula= as.formula(paste(colnames(mydat@data)[i+2], "~", "cov1 + cov2")), 
    coordinates=mydat@coords,
    params=paramssn, 
    paramToEstimate = c('nugget','anisoRatio'),
    boxcox = 1,
    cilevel=0.8,
    type = "double",
    reml=FALSE, 
    NparamPerIter=400,
    Nglobal=c(128,64),
    Nlocal=c(16,16),
    NlocalCache=2800,
    verbose=FALSE)
  
  a <- geostatsp::loglikLgm(
    c(trueParam['range'],
      trueParam['shape'],
      nugget,
      trueParam[c('anisoRatio')],
      trueParam[c('anisoAngleRadians')],
      boxcox),
    data = mydat,
    formula = as.formula(paste(colnames(mydat@data)[i+2], "~", "cov1 + cov2")),
    reml = FALSE,
    minustwotimes=FALSE)[['logLik']]
  
  if (a  >= result$breaks){
    count = count +1
  }
  breaksf[i] = result$breaks
  makedata0 <- cbind(paramssn[,1:2], result$LogLik, D=i)
  makedata <- rbind(makedata, makedata0)
  estimates <- cbind(estimates, result$summary)
}
colnames(makedata) <- c('nugget','anisoRatio', 'LogLik', 'D')
makedatahey <- makedata[-1,]
estimates_summary <- estimates[,-1]
mean_estimates <- apply(estimates_summary, 1, mean)
count/15
mean_estimates




library(RColorBrewer)
library('ggplot2')
#rm(forigin_ns)
forigin_ns <- ggplot(subset(makedatahey, D==1), aes(x = nugget, y = anisoRatio, z=LogLik)) + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_point(size = 0.5, colour='grey') + 
  stat_contour(breaks=breaksf[1])  + 
  geom_vline(xintercept = 4)+
  geom_hline(yintercept=trueParam['anisoRatio'])

cols <- rainbow(15)
for (i in 2:15){
  forigin_ns <- forigin_ns + geom_contour(data = subset(makedatahey, D==i), aes(x = nugget, y = anisoRatio, z=LogLik), breaks=breaksf[i], col= cols[i])
}
forigin_ns


# > count/15
# [1] 0.8666667
# > mean_estimates
# (Intercept)        cov1        cov2   sdSpatial      nugget  anisoRatio 
# 3.0759438   0.9996654   0.5005581   0.5319227   3.9809522   4.3333333 







