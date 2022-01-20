############### A SIMULATION STUDY IN THE PAPER 2 ####################################
library('geostatsp')
Nsim=10
Ngrid=20
set.seed(88)
# pointsGrid = expand.grid(seq(100*1e4,120*1e4, len=Ngrid), seq(100*1e4,120*1e4, len=Ngrid))
# #plot(pointsGrid,cex=0.3)
# pointsGrid = pointsGrid[,1] + 1i*pointsGrid[,2]
# pointsRandom = pointsGrid + stats::runif(length(pointsGrid), 500, 600) * exp(stats::runif(length(pointsGrid), 0, 2*pi)*1i)
# plot(pointsRandom,cex=0.3)
# #plot(Re(pointsRandom),Im(pointsRandom),cex=0.3)
# mydat = SpatialPointsDataFrame(cbind(Re(pointsRandom),Im(pointsRandom)),
#                                data=data.frame(cov1 =stats::rnorm(Ngrid^2, 120, 20), cov2 = stats::rpois(Ngrid^2, 50))
# )


mydat = SpatialPointsDataFrame(cbind(seq(100*1e3,140*1e3, len=500), seq(100*1e3,120*1e3, len=500)),
                               data=data.frame(cov1 =stats::rnorm(500, 120, 20), cov2 = stats::rpois(500, 50))
)
plot(cbind(seq(100*1e3,140*1e3, len=500), seq(100*1e3,120*1e3, len=500)), xlab="",ylab="")

# par(mfrow = c(1, 2))
# library('geoR')
# sim1 <- grf(400, cov.pars = c(1, 0.25))
# points.geodata(sim1, main = "simulated locations and values")
# plot(sim1, max.dist = 1, main = "true and empirical variograms")


#mydat@data
#mydat@coords

## simulate a random field
trueParam = c(range=1000, shape=2, variance=0.5^2, tauSq=1^2, anisoRatio=2, anisoAngleRadians=0.3)
trueParam['nugget'] = trueParam['tauSq']/trueParam['variance']
boxcox = c(boxcox=1)
options("useRandomFields" = FALSE)
for (i in 1:Nsim){
  set.seed(i*6-1)
  U <- geostatsp::RFsimulate(model=trueParam[setdiff(names(trueParam),'tauSq')],x=mydat)@data
  # simulate response
  Y <- 3 + 1*mydat$cov1 + 0.5*mydat$cov2 + U + stats::rnorm(length(mydat), 0, sd=sqrt(trueParam["tauSq"]))
  mydat@data <- cbind(mydat@data,Y)
}
colnames(mydat@data) <- c("cov1", "cov2", paste("Y", c(1:Nsim), sep=""))
#mydat@data[1:20,]
#mydat$Y3

## geostatsp's estimates
swissRes =  lgm( formula=Y1~ cov1 + cov2,
                 data=mydat,
                 grid=20,
                 reml = FALSE,
                 shape=2,
                 boxcox=1,
                 fixBoxcox=TRUE, fixShape=TRUE, fixNugget = FALSE,  #Set to FALSE to estimate the nugget effect parameter.
                 aniso=TRUE )
swissRes$summary[,c('estimate','ci0.05', 'ci0.95')]
swissRes$optim$mle



## ---gpuRandom------range vs nugget contour--------------------------------------------------------------
## set params
newParamList = list(
  range = c(exp(seq(log(0.2), log(1), len=10))*1000, exp(seq(log(1.2), log(4), len=10))*1000),  #22
  nugget = c(seq(0, 4, len=8), seq(4.5, 20, len=12)),  #20
  shape = 2,
  anisoRatio =2,
  anisoAngleRadians = 0.3
) 
params = do.call(expand.grid, newParamList)

library(gpuRandom)
library(gpuR)
count = 0
breaksf = rep(0, 10)
makedata <- rep(0, 4)
estimates <- rep(0, 6)
## calculate logl
for (i in 1:5){
  result<-gpuRandom::likfitLgmCov2d(
    data= mydat,
    formula= as.formula(paste(colnames(mydat@data)[i+2], "~", "cov1 + cov2")), 
    coordinates=mydat@coords,
    params=params, 
    paramToEstimate = c('range','nugget'),
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
      trueParam['nugget'],
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
  makedata0 <- cbind(params[,1:2], result$LogLik, D=i)
  makedata <- rbind(makedata, makedata0)
  estimates <- cbind(estimates, result$summary)
}


colnames(makedata) <- c("range",'nugget', 'LogLik', 'D')
makedatahey <- makedata[-1,]
estimates_summary <- estimates[,-1]
mean_estimates <- apply(estimates_summary, 1, mean)
count/5
mean_estimates



library(RColorBrewer)
library('ggplot2')
#rm(forigin)
forigin <- ggplot(subset(makedatahey, D==1), aes(x = range, y = nugget, z=LogLik)) + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_point(size = 0.5, colour='grey') + 
  stat_contour(breaks=breaksf[1])  + 
  geom_vline(xintercept = trueParam['range'])+
  geom_hline(yintercept=trueParam['nugget'])


cols <- rainbow(5)
for (i in 2:5){
  forigin <- forigin + geom_contour(data = subset(makedatahey, D==i), aes(x = range, y = nugget, z=LogLik), breaks=breaksf[i], col= cols[i])
}
forigin
# (Intercept)         cov1         cov2    sdSpatial        range       nugget 
# 3.0807703    0.9996489    0.5005081    0.5338107 1040.7750111    3.9151560 
# 
# i=1
a <- subset(makedatahey, D==i)[,3]
lMatrix = matrix(a, length(newParamList[[1]]), length(newParamList[[2]]))
contour(newParamList[[1]], newParamList[[2]], lMatrix,
        col = par("fg"), lty = par("lty"), lwd = par("lwd"),
        add = FALSE, levels = c(
                                breaksf[i]),  xlab = "range", ylab = "nugget")

abline(h=2, col='red')
abline(v=1000, col='red')
points(x = params[,1], y = params[,2], pch=20, col='grey', cex=0.5)









## ---gpuRandom------range vs shape contour--------------------------------------------------------------
## set params
newParamListrs = list(
    range = c(exp(seq(log(0.2), log(1), len=8))*1000, exp(seq(log(1.2), log(4), len=12))*1000),  
    shape = c(0.1, 0.15, 0.2, 0.23, 0.3, 0.4, 0.5, 0.75, 0.9, 1, 2, 2.5, 3, 4),
    nugget = 4,  
    anisoRatio =1,
    anisoAngleRadians = 0.0
  ) 
paramsrs = do.call(expand.grid, newParamListrs)

count = 0
breaksf = rep(0, 10)
makedata <- rep(0, 4)
estimates <- rep(0, 6)
## calculate logl
for (i in 1:10){
  result<-gpuRandom::likfitLgmCov2d(
    data= mydat,
    formula= as.formula(paste(colnames(mydat@data)[i+2], "~", "cov1 + cov2")), 
    coordinates=mydat@coords,
    params=paramsrs, 
    paramToEstimate = c('range','shape'),
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
  makedata0 <- cbind(paramsrs[,1:2], result$LogLik, D=i)
  makedata <- rbind(makedata, makedata0)
  estimates <- cbind(estimates, result$summary)
}
colnames(makedata) <- c("range",'shape', 'LogLik', 'D')
makedatahey <- makedata[-1,]
estimates_summary <- estimates[,-1]
mean_estimates <- apply(estimates_summary, 1, mean)
count/15
mean_estimates


library(RColorBrewer)
library('ggplot2')
#rm(forigin_s)
forigin_s <- ggplot(subset(makedatahey, D==1), aes(x = range, y = shape, z=LogLik)) + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_point(size = 0.5, colour='grey') + 
  stat_contour(breaks=breaksf[1])  + 
  geom_vline(xintercept = trueParam['range'])+
  geom_hline(yintercept=trueParam['shape'])

cols <- rainbow(15)
for (i in 2:15){
  forigin_s <- forigin_s + geom_contour(data = subset(makedatahey, D==i), aes(x = range, y = shape, z=LogLik), breaks=breaksf[i], col= cols[i])
}
forigin_s


par(mfrow = c(2, 3))
i=2
`a <- subset(makedatahey, D==i)[,3]
lMatrix = matrix(a, length(newParamListrs[[1]]), length(newParamListrs[[2]]))
contour(newParamListrs[[1]], newParamListrs[[2]], lMatrix,
        col = par("fg"), lty = par("lty"), lwd = par("lwd"),
        add = FALSE, levels = c(breaksf[i]-4, breaksf[i]-3, 
                                breaksf[i]-2, breaksf[i]-1,
                                breaksf[i]),  xlab = "range", ylab = "shape")

abline(h=2, col='red')
abline(v=1000, col='red')
points(x = paramsrs[,1], y = paramsrs[,2], pch=20, col='grey', cex=0.5)
`
####################################################################







########################################################################################################################################
## ---gpuRandom------nugget vs shape contour--------------------------------------------------------------
## set params
newParamListsn = list(
  nugget = c(seq(0,4,len=8), seq(5,20, len=8)),  
  shape = c(0.1, 0.15, 0.2, 0.23, 0.3, 0.4, 0.45, 0.5, 0.75, 0.9, 1, 2, 3, 5),
  range = 1000,  
  anisoRatio =1,
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
    paramToEstimate = c('nugget','shape'),
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
colnames(makedata) <- c("nugget",'shape', 'LogLik', 'D')
makedatahey <- makedata[-1,]
estimates_summary <- estimates[,-1]
mean_estimates <- apply(estimates_summary, 1, mean)
count/15
mean_estimates


library(RColorBrewer)
library('ggplot2')
#rm(forigin_ns)
forigin_ns <- ggplot(subset(makedatahey, D==1), aes(x = nugget, y = shape, z=LogLik)) + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_point(size = 0.5, colour='grey') + 
  stat_contour(breaks=breaksf[1])  + 
  geom_vline(xintercept = nugget)+
  geom_hline(yintercept=trueParam['shape'])

cols <- rainbow(15)
for (i in 2:15){
  forigin_ns <- forigin_ns + geom_contour(data = subset(makedatahey, D==i), aes(x = nugget, y = shape, z=LogLik), breaks=breaksf[i], col= cols[i])
}
forigin_ns









########################################################################################################################################
## ---gpuRandom------nugget vs shape contour--------------------------------------------------------------
## set params
newParamListsn = list(
  nugget = c(seq(0,4,len=8), seq(5,20, len=8)),  
  shape = c(0.1, 0.15, 0.2, 0.23, 0.3, 0.4, 0.45, 0.5, 0.75, 0.9, 1, 2, 3, 5),
  range = 1000,  
  anisoRatio =1,
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
    paramToEstimate = c('nugget','shape'),
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
colnames(makedata) <- c("nugget",'shape', 'LogLik', 'D')
makedatahey <- makedata[-1,]
estimates_summary <- estimates[,-1]
mean_estimates <- apply(estimates_summary, 1, mean)
count/15
mean_estimates


library(RColorBrewer)
library('ggplot2')
#rm(forigin_ns)
forigin_ns <- ggplot(subset(makedatahey, D==1), aes(x = nugget, y = shape, z=LogLik)) + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_point(size = 0.5, colour='grey') + 
  stat_contour(breaks=breaksf[1])  + 
  geom_vline(xintercept = nugget)+
  geom_hline(yintercept=trueParam['shape'])

cols <- rainbow(15)
for (i in 2:15){
  forigin_ns <- forigin_ns + geom_contour(data = subset(makedatahey, D==i), aes(x = nugget, y = shape, z=LogLik), breaks=breaksf[i], col= cols[i])
}
forigin_ns














