########################################################################################################################################
## ---gpuRandom------shape vs ratio contour--------------------------------------------------------------
## set params
newParamListsn = list(
  shape = c(0.1, 0.15, 0.2, 0.23, 0.3, 0.4, 0.45, 0.5, 0.75, 0.9, 1, 2, 3, 5, 6),
  anisoRatio =c(seq(0.1, 1, len=2), seq(2,10, len=8)),
  nugget = 4,  
  range = 1000,  
  anisoAngleRadians = 0.0
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
    paramToEstimate = c('shape','anisoRatio'),
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
colnames(makedata) <- c('shape',"anisoRatio", 'LogLik', 'D')
makedatahey <- makedata[-1,]
estimates_summary <- estimates[,-1]
mean_estimates <- apply(estimates_summary, 1, mean)
count/15
mean_estimates


library(RColorBrewer)
library('ggplot2')
#rm(forigin_ns)
forigin_ns <- ggplot(subset(makedatahey, D==1), aes(x = shape, y = anisoRatio, z=LogLik)) + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_point(size = 0.5, colour='grey') + 
  stat_contour(breaks=breaksf[1])  + 
  geom_vline(xintercept = trueParam['shape'])+
  geom_hline(yintercept=trueParam['anisoRatio'])

cols <- rainbow(15)
for (i in 2:15){
  forigin_ns <- forigin_ns + geom_contour(data = subset(makedatahey, D==i), aes(x = shape, y = anisoRatio, z=LogLik), breaks=breaksf[i], col= cols[i])
}
forigin_ns












