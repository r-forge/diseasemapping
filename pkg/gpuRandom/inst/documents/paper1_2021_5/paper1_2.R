## check if needee packages are installed
packagesList <- c("knitr", "tidyverse", "kableExtra", "devtools")
for (i in 1:length(packagesList)){
  if (!requireNamespace(packagesList[i], quietly = TRUE))
    install.packages(packagesList[i])
}

if (!requireNamespace("gpuR", quietly = TRUE))
  devtools::install_github("cdeterman/gpuR")

if (!requireNamespace("geostatsp", quietly = TRUE))
  install.packages("geostatsp", repos='http://r-forge.r-project.org')

if (!requireNamespace("clrng", quietly = TRUE))
  devtools::install_github("ruoyongxu/clrng")

if (!requireNamespace("gpuBatchMatrix", quietly = TRUE))
  devtools::install_github("ruoyongxu/gpuBatchMatrix")




library("gpuR")
setContext(   grep('gpu', listContexts()$device_type) [1]    )

sessionInfo()
gpuInfo()

## Section 2.1
# Creating streams on CPU
library("clrng")
setBaseCreator(rep(12345,6))
myStreamsCpu <- createStreamsCpu(4)
t(myStreamsCpu)

# Creating streams on GPU
myStreamsGpu = vclMatrix(myStreamsCpu)

## Section 2.2
# Generate 6 i.i.d. U (0,1) random numbers
sim_1 =  runifGpu(n = 8, streams = myStreamsGpu, Nglobal = c(2,2))
as.vector(sim_1)

## Section 3.1
# Generate a large matrix of normal random numbers, test the run time
streams <- createStreamsGpu(512 * 128)
system.time(rnormGpu(c(10000,10000), streams=streams, Nglobal=c(512,128), type="double"))
system.time(matrix(rnorm(10000^2),10000,10000))


## Section 3.2
# Generate exponential random numbers
r_matrix <- rexpGpu(c(2,4), rate=1, myStreamsGpu, Nglobal=c(2,2), type="double")
as.matrix(r_matrix)


## Section 4.1
data("month")
month <- as.matrix(month)
library(knitr)
kable(month, format = "markdown", caption = "Monthly birth anomaly data")




# using GPU
streams <- createStreamsGpu(n = 256*64)
month_gpu <- vclMatrix(month, type = "integer")
system.time(result_month <- clrng::fisher.sim(month_gpu, 1e6, streams=streams, type="double", returnStatistics=TRUE,  Nglobal = c(256,64)))

unlist(result_month[c('threshold','simNum','counts')])
result_month$p.value


# save test statistics for plot
month_stats <- as.vector(result_month$sim)

# using CPU
# notice this can take about more than 1 minute
system.time(result_monthcpu<-stats::fisher.test(month,simulate.p.value = TRUE, B=1015808))
result_monthcpu$p.value


## Section 4.2
data("week")
week <- as.matrix(week)
kable(week, format = "markdown", caption = "Day-of-week birth anomaly data") 

# using GPU
week_GPU<-gpuR::vclMatrix(week,type="integer")
system.time(resultWeek<-clrng::fisher.sim(week_GPU, 1e7, streams=streams, type="double",returnStatistics=TRUE,Nglobal = c(256,64)))
unlist(resultWeek[c('threshold','simNum','counts')])
resultWeek$p.value


# save test statistics for plot
resultWeek$cpu = as.vector(resultWeek$sim)


# using CPU
# notice that running in rstudio can take about 10 minutes 
# set a smaller number for B if you think that's too long time
system.time(result_weekcpu<-fisher.test(week,simulate.p.value = TRUE,B=10010624))
result_weekcpu$p.value


## summary table
library(kableExtra)
dt <- data.frame(
  var1 = c("B", "1M", "10M", "1M", "10M"),
  var2 = c('Intel 2.5ghz', 0.403804, 0.0001251, 10.74,  63.24),
  var3 = c('AMD Radeon', 0.403507, 0.0001274,   2.28,    10.82),
  var4 = c('Intel 3.7ghz', 0.4035606,0.0001202, 15.00,  91.58),
  var5 = c('NVIDIA V100', 0.403507,0.0001274,   0.72,   4.18),
  var6 = c('Data','month','week','month','week')
)

knitr::kable(dt, col.names = NULL, caption = "Summary of comparions of Fisher's test simulation on different devices. Computer 1 is equipped with CPU Intel Xenon W-2145 3.7Ghz and AMD Radeon VII. 
             Computer 2 is equipped with VCPU Intel Xenon Skylake 2.5Ghz and VGPU Nvidia Tesla V100.") %>%
  kable_styling(full_width = FALSE, position = "center") %>%
  kableExtra::group_rows(index = c("", "P-value" = 2, "Run-time" = 2)) %>%
  add_header_above(c(" ", "Computer 1" = 2, "Computer 2" = 2, " " = 1))


## test statistics plot
hist(month_stats, xlab="test statistic",#TeX('-$\\sum(\\log(n_{ij}!))$'),
     breaks=40, ylab="proportion", main="", prob=TRUE)
abline(v = result_month$threshold, col = "blue", lwd = 1.5)

hist(resultWeek$cpu, xlab='test statistic', ylab="proportion", prob=TRUE, main="")
abline(v = resultWeek$threshold, col = "blue", lwd = 1.5)


## Section 5.1
library("gpuBatchMatrix")
library('geostatsp')


# setupcoords
NlocalCache = 1000
Nglobal = c(128, 64, 2)
Nlocal = c(16, 8, 2)
theType = "double"


# swissRainBoundary
data("swissRain", package="geostatsp")
myRaster = geostatsp::squareRaster(swissBorder, 90)
dim(myRaster)


# setup paramsBatch
params = 
  rbind(c(shape=1.25, range=50*1000, variance = 1.5, anisoRatio = 1, anisoAngleRadians = 0), 
        c(2.15, 60*1000, 2,  4, pi/7),
        c(0.6, 30*1000, 2,  2, pi/5),
        c(3, 30*1000, 2,   2, pi/7)
  )
params



# matern batches
maternCov = gpuBatchMatrix::maternBatch(
  params, myRaster,          
  Nglobal=c(128,64), Nlocal=c(16,4))
dim(maternCov)


# take cholesky decomposition
diagMat = gpuBatchMatrix::cholBatch(maternCov, Nglobal= c(128, 8), Nlocal= c(32, 8))


# simulate random normal vectors
streamsForGrf <- createStreamsGpu(128*64)
zmatGpu = clrng::rnormGpu(
                c(nrow(maternCov),2), streams=streamsForGrf, 
                Nglobal=c(128,64),
                type = 'double')


# L*D*Z
simMat = gpuBatchMatrix::multiplyLowerDiagonalBatch(maternCov, 
                                                    diagMat, zmatGpu,
                                                    diagIsOne = TRUE,   
                                                    transformD = "sqrt", 
                                                    Nglobal=c(128, 64, 2), 
                                                    Nlocal= c(8, 2, 1), 
                                                    NlocalCache=1000)


# plot Setup
simRaster = raster::brick(myRaster, nl = ncol(simMat)*nrow(params))
values(simRaster) = as.vector(as.matrix(simMat))
plot(simRaster)
names(simRaster) = apply(expand.grid('par',1:nrow(params), 
                                     'sim', 1:ncol(simMat)), 1, paste, collapse='')
theSubcap = gsub("par", "parameter ", names(simRaster))
theSubcap = gsub("sim", ", simuation ", theSubcap)

# random surfaces plot
myCol = mapmisc::colourScale(breaks=c(-7, -4, seq(-2, 2), 4, 7), style='fixed', col='Spectral')
for(D in names(simRaster)) {
  mapmisc::map.new(simRaster)
  plot(simRaster[[D]], legend=FALSE, add=TRUE, col=myCol$col, breaks=myCol$breaks)
  plot(swissBorder, add=TRUE)
}
mapmisc::legendBreaks("topleft", myCol, inset=0)


sessionInfo()
gpuInfo()









