############### A SIMULATION STUDY IN THE PAPER 2 ####################################
library('geostatsp')
Nsim=10
set.seed(88)
mydat = SpatialPointsDataFrame(cbind(seq(100*1e3,140*1e3, len=500), seq(100*1e3,120*1e3, len=500)), 
                               data=data.frame(cov1 =stats::rnorm(500, 120, 20), cov2 = stats::rpois(500, 50))
)
#mydat@data
#mydat@coords

## simulate a random field
trueParam = c(range=1000, shape=2, variance=0.5^2, tauSq=1^2, anisoRatio=1, anisoAngleRadians=0)
nugget = trueParam['tauSq']/trueParam['variance']
names(nugget) <- "nugget"
nugget   #4
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
swissRes =  lgm( formula=Y10~ cov1 + cov2,
                 data=mydat,
                 grid=20,
                 reml = FALSE,
                 boxcox=1,
                 fixBoxcox=TRUE, fixShape=FALSE, fixNugget = FALSE,  #Set to FALSE to estimate the nugget effect parameter.
                 aniso=FALSE )
swissRes$summary[,c('estimate','ci0.005', 'ci0.995')]
swissRes$optim$mle





## ---gpuRandom------range vs nugget conntour--------------------------------------------------------------
## set params
newParamList = list(
  range = c(exp(seq(log(0.3), log(1), len=10))*1000, exp(seq(log(1.2), log(2.9), len=12))*1000),  #22
  nugget = c(seq(0, 4, len=8), seq(4.5, 20, len=12)),  #20
  shape = 2,
  anisoRatio =1,
  anisoAngleRadians = 0.0
) 
params = do.call(expand.grid, newParamList)

result<-gpuRandom::likfitLgmCov2d(
  data= mydat,
  formula= as.formula(paste(colnames(mydat@data)[i+2], "~", "cov1 + cov2")), 
  coordinates=mydat@coords,
  params=params, # CPU matrix for now, users need to provide proper parameters given their specific need
  paramToEstimate = c('range','nugget'),
  boxcox = c(1,0),# (seq(-1,3,len=8),
  cilevel=0.8,
  type = "double",
  reml=FALSE, 
  NparamPerIter=400,
  Nglobal=c(128,64),
  Nlocal=c(16,16),
  NlocalCache=2800,
  verbose=FALSE)

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
    params=params, # CPU matrix for now, users need to provide proper parameters given their specific need
    paramToEstimate = c('range','nugget'),
    boxcox = 1,# (seq(-1,3,len=8),
    cilevel=0.8,
    type = "double",
    reml=FALSE, 
    NparamPerIter=400,
    Nglobal=c(128,64),
    Nlocal=c(16,16),
    NlocalCache=2800,
    verbose=FALSE)

  breaksf[i] = result$breaks
  LogLik = result$LogLik[,1]
  #max(result$LogLik[,1])-qchisq(0.8,  df = 2)/2
  #as.vector(result$jacobian)
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
  makedata0 <- cbind(params[,1:2], result$LogLik[,1], D=i)
  makedata <- rbind(makedata, makedata0)
  estimates <- cbind(estimates, result$summary)
}




count/Nsim
colnames(makedata) <- c("range",'nugget', 'LogLik', 'D')
makedatahey <- makedata[-1,]
estimates_summary <- estimates[,-1]

library(RColorBrewer)
library('ggplot2')
#rm(forigin)
forigin <- ggplot(subset(makedatahey, D==1), aes(x = range, y = nugget, z=LogLik)) + theme_bw() + 
  scale_y_continuous(limits=c(0,20))+
  scale_x_continuous(limits=c(300,3000) )+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_point(size = 0.5, colour='grey') + 
  stat_contour(breaks=breaksf[1])  + 
  geom_vline(xintercept = trueParam['range'])+
  geom_hline(yintercept=nugget)

# cols <- c("#FFDB6D", "#C4961A", "#F4EDCA", "#CC79A7","#999999",
#           "#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352")
# cols <- function(n) hcl.colors(n, "ag_Sunset")
# #ggsave("filename", plot = myPlot)
# cols <- cols(Nsim)
# cols <- brewer.pal(14,'Set5')
cols <- rainbow(Nsim/2)
for (i in 2:Nsim/2){
  forigin <- forigin + geom_contour(data = subset(makedatahey, D==i), aes(x = range, y = nugget, z=LogLik), breaks=breaksf[i], col= cols[i])
}
forigin











