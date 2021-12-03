## check if needee packages are installed
packagesList <- c("knitr", "tidyverse", "kableExtra", "devtools")

for (i in 1:length(packagesList)){
  if(packagesList[i] %in% rownames(installed.packages()) == FALSE) 
   {install.packages(packagesList[i])}
}

if("gpuR" %in% rownames(installed.packages()) == FALSE) 
{devtools::install_github("cdeterman/gpuR")}    #{install.packages("gpuR", repos='https://github.com/cdeterman/gpuR')}

if("geostatsp" %in% rownames(installed.packages()) == FALSE) 
{install.packages("geostatsp", repos='http://r-forge.r-project.org')}

if("clrng" %in% rownames(installed.packages()) == FALSE) 
{devtools::install_github("Ruoyong/clrng")}

if("gpuBatchMatrix" %in% rownames(installed.packages()) == FALSE) 
{devtools::install_github("Ruoyong/gpuBatchMatrix")}




## ----downloadJssClassFile, eval=FALSE, include=FALSE-------------------------------------------------------------------------------------------------------
## zUrl = 'http://www.jstatsoft.org/public/journals/1/jss-style.zip'
## zFile = basename(zUrl)
## if(!file.exists(zFile))
##         download.file(zUrl, zFile)
## unzip(zFile, exdir='.',
##                 files = grep("^jss", unzip(zFile, list=TRUE)$Name, value=TRUE)
## )


## ----echo=FALSE,results="hide"-----------------------------------------------------------------------------------------------------------------------------
options(width=65)


## ----setup, cache=FALSE,include=FALSE----------------------------------------------------------------------------------------------------------------------
library('knitr')
options(continue="+  ",  prompt="R> ", digits=3, width=57) 

opts_chunk$set(echo=TRUE,fig.path='Figures/G', fig.height=3,
                fig.width=4.5,dev='png',tidy=TRUE,fig.align='center',
                tidy.opts=list(blank=FALSE, width.cutoff=57),
                prompt=TRUE, highlight=FALSE,cache.stuff=1)

hook_output <- function(x, options) {
        if (knitr:::output_asis(x, options)) return(x)
        paste0('\\begin{CodeOutput}\n', x, '\\end{CodeOutput}\n')
}

knit_hooks$set(output = hook_output)
knit_hooks$set(message = hook_output)
knit_hooks$set(warning = hook_output)


knit_hooks$set(source  =function(x, options) {
        x = gsub("\n", "\n+", x)
        x = paste("R>",paste(x, collapse='\nR> '))
#       paste0(c('\\begin{example*}', x, '\\end{example*}', ''),
#                       collapse = '\n')
x =  paste0(c('\\begin{CodeInput}', x, '\\end{CodeInput}', ''),
                        collapse = '\n')
                        x               
})

knit_hooks$set(plot = function (x, options) 
                {
                        paste( hook_plot_tex(x, options), "\n", 
                                        sep = "")
                }
)

knit_hooks$set(chunk  = function(x, options) {
        if (knitr:::output_asis(x, options)) return(x)
        theend = gregexpr("end\\{Code(Input|Output)\\}", x)
        theend=theend[[1]]
        if(theend[1]>0){
                theend = theend[length(theend)] + attributes(theend)$match.length[length(theend)]
                x = paste(substr(x,1,theend-1),"\n\\end{CodeChunk}", substr(x,theend, nchar(x)))
                x=sub('begin\\{Code', 'begin{CodeChunk}\n\\\\begin{Code', x)
        }
        x
})
knitr::knit_hooks$set(margins = function(before, options, envir) {
    if (!before) 
        return()
    graphics::par(mar = c(1.5 + 0.9 * options$margins, 1.5 + 
        0.9 * options$margins, 0.2, 0.2), mgp = c(1.45, 0.45, 
        0), cex = 1.25)
})
opts_knit$set(header='')



## ----preliminaries, echo=FALSE, results="hide", eval=FALSE-------------------------------------------------------------------------------------------------
##   knitr::knit_hooks$set(plot=knitr::hook_plot_tex)
##   options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)
## knitr::opts_chunk$set(highlight=FALSE, background='#FFFFFF00', fig.height=4, fig.width=6, out.width='0.45\\textwidth')
## library("MASS")


## ----echo=FALSE, results="hide", message=FALSE-------------------------------------------------------------------------------------------------------------
library("gpuR")
library("clrng")
setContext(   grep('gpu', listContexts()$device_type) [1]    )
#setwd("/home/ruoyong/diseasemapping/pkg/clrng/inst/documents/paper1_2021_5")


## ----createStreams_CPU-------------------------------------------------------------------------------------------------------------------------------------
# creating streams on CPU
library("gpuR")
library("clrng")
myStreamsCpu <- createStreamsCpu(n=4, initial=12345)
t(myStreamsCpu)


## ----createStreams_Gpu-------------------------------------------------------------------------------------------------------------------------------------
# creating streams on GPU
myStreamsGpu = vclMatrix(myStreamsCpu)
myStreamsGpu2 = createStreamsGpu(n=4, initial=12345)


## ----simStuff2, eval=TRUE,tidy=TRUE, echo=TRUE-------------------------------------------------------------------------------------------------------------
as.vector(clrng::runif(n=6, streams=myStreamsGpu, Nglobal=c(2,2)))


## ----eval=TRUE,tidy=TRUE, include=FALSE--------------------------------------------------------------------------------------------------------------------
t(matrix(as.matrix(myStreamsGpu), nrow(myStreamsCpu), ncol(myStreamsCpu), dimnames = dimnames(myStreamsCpu)))


## ----save streams,echo=TRUE,eval=TRUE----------------------------------------------------------------------------------------------------------------------
saveRDS(as.matrix(createStreamsGpu(n=4)), "myStreams.rds")
# Load the streams object as streams_saved
streams_saved <- vclMatrix(readRDS("myStreams.rds"))


## ----random normal time compare, eval=TRUE, cache=FALSE----------------------------------------------------------------------------------------------------
streams <- createStreamsGpu(n = 512 * 128)
system.time(clrng::rnorm(c(10000,10000), streams=streams, 
                             Nglobal=c(512,128), type="double"))

## ----random normal time compare2, eval=TRUE, cache=TRUE----------------------------------------------------------------------------------------------------
system.time(matrix(stats::rnorm(10000^2),10000,10000))


## ----random exponential, eval=TRUE, cache=FALSE------------------------------------------------------------------------------------------------------------
r_matrix <- clrng::rexp(c(2,4), rate=1, Nglobal=c(2,2), type="double")
as.matrix(r_matrix)


## ----Timecomparemonth,echo=FALSE,eval=TRUE, cache=FALSE, results="hide"------------------------------------------------------------------------------------
datamonth<-read.csv("/home/ruoyong/diseasemapping/pkg/gpuRandom/inst/documents/paper1_2021_5/month.csv") 
month<-as.matrix(datamonth[,-1])
rownames(month) <- c("Jan", "Feb","Mar","Apr","May","Jun","Jul","Aug","Sept","Oct","Nov","Dec")
colnames(month) <- c("Ane", "Men", "Cya", "Her", "Omp", "Gas", "Lim", "Cle", "Pal", "Dow", "Chr", "Hyp")


## ----monthdata,eval=TRUE,echo=FALSE, cache=TRUE, fig.pos='h', message=FALSE--------------------------------------------------------------------------------
library(knitr)
library(tidyverse)
library(kableExtra)
knitr::kable(month,format="latex", align = c("rrrrrrrrrrrr"), #label="tab:month"
 caption = "Monthly birth anomaly data\\label{tab:month}", booktabs=TRUE,linesep = "") %>%
 kable_styling(full_width = F, position = "center")#latex_options = "HOLD_position")


## ----timecomparemonth1Gpu, eval=TRUE,cache=FALSE-----------------------------------------------------------------------------------------------------------
streams <- createStreamsGpu(n =256*64, initial=666)
month_gpu<-vclMatrix(month,type="integer")
system.time(result_month <- clrng::fisher.sim(month_gpu, 1e6, streams=streams,
                   type="double", returnStatistics=TRUE,  Nglobal = c(256,64)))
result_month$threshold
result_month$simNum
result_month$counts
result_month$p.value


## ----timecomparemonth3Gpu, eval=TRUE,echo=FALSE, cache=TRUE------------------------------------------------------------------------------------------------
month_stats <- as.vector(result_month$sim)


## ----timecomparemonth4Cpu, cache=TRUE, eval=TRUE-----------------------------------------------------------------------------------------------------------
#using CPU
system.time(result_monthcpu<-stats::fisher.test(month,simulate.p.value = TRUE, B=1e6))
result_monthcpu$p.value


## ----Timecompareweek, echo=FALSE, eval=TRUE, cache=TRUE,results="hide", fig.pos='h'------------------------------------------------------------------------
dataweek<-read.csv("/home/ruoyong/diseasemapping/pkg/gpuRandom/inst/documents/paper1_2021_5/weekday.csv") 
week<-as.matrix(dataweek[,-1])
rownames(week) <- c("Mon", "Tue","Wed","Thu","Fri","Sat","Sun")
colnames(week) <- c("Ane", "Men", "Cya", "Her", "Omp", "Gas", "Lim", "Cle", "Pal", "Dow", "Chr", "Hyp")


## ----weekdata, eval=TRUE, echo=FALSE, cache=TRUE-----------------------------------------------------------------------------------------------------------
knitr::kable(week,format="latex",align = c("rrrrrrrrrrrr"),  #label = "tab:week",
             caption = "Day-of-week birth anomaly data\\label{tab:week}", booktabs=TRUE,linesep = "") %>%
kable_styling(full_width = F, position = "center")#, latex_options = "HOLD_position")


## ----TimecompareweekGpu2, eval=TRUE, cache=FALSE-----------------------------------------------------------------------------------------------------------
week_GPU<-gpuR::vclMatrix(week,type="integer")
system.time(result_week<-clrng::fisher.sim(week_GPU, 1e7, streams=streams,
                type="double",returnStatistics=TRUE,Nglobal = c(256,64)))
result_week$threshold
result_week$simNum
result_week$counts
result_week$p.value


## ----TimecompareweekGpu3, echo=FALSE, eval=TRUE, cache=TRUE------------------------------------------------------------------------------------------------
result_week$cpu = as.vector(result_week$sim)


## ----TimecompareweekCpu, cache=TRUE, eval=TRUE-------------------------------------------------------------------------------------------------------------
#using CPU
system.time(result_weekcpu<-fisher.test(week,simulate.p.value = TRUE,B=10010624))
result_weekcpu$p.value


## ----summarycompare, echo=FALSE, eval=TRUE, cache=TRUE, fig.pos='h'----------------------------------------------------------------------------------------
library(kableExtra)
dt <- data.frame(
  var1 = c("B", "1M", "10M", "1M", "10M"),
  var2 = c('Intel 2.5ghz', 0.403804, 0.0001251, 10.74,  63.24),
  var3 = c('AMD Radeon', 0.403507, 0.0001274,   2.28,    10.82),
  var4 = c('Intel 3.7ghz', 0.4035606,0.0001202, 15.00,  91.58),
  var5 = c('NVIDIA V100', 0.403507,0.0001274,   0.72,   4.18),
  var6 = c('Data','month','week','month','week')
)

knitr::kable(dt, col.names = NULL, caption = "Summary of comparions of Fisher's test simulation on different devices. Computer 1 is equipped with CPU Intel Xenon W-2145 3.7Ghz and AMD Radeon VII. Computer 2 is equipped with VCPU Intel Xenon Skylake 2.5Ghz and VGPU Nvidia Tesla V100.\\label{tab:summary}") %>%
  kable_styling(full_width = FALSE, position = "center") %>%
  kableExtra::group_rows(index = c("P-value" = 2, "Run-time" = 2)) %>%
  add_header_above(c(" " = 1, "Computer 1" = 2, "Computer 2" = 2, " " = 1))


## ----fighistMonth, eval=TRUE, echo=FALSE, dev='pdf', fig.cap = "Approximate sampling distributions of the test statistics from the two examples Month and Week. The test statistic values of the observed tables are marked with a blue line on each plot.", fig.subcap = c('Month data', 'week data'), out.width="0.47\\textwidth", margins=1----
hist(month_stats, xlab="test statistic",#TeX('-$\\sum(\\log(n_{ij}!))$'),
     breaks=40, ylab="proportion", main="", prob=TRUE)
abline(v = result_month$threshold, col = "blue", lwd = 1.5)

hist(result_week$cpu, xlab='test statistic', ylab="proportion", prob=TRUE, main="")
abline(v = result_week$threshold, col = "blue", lwd = 1.5)


## ----setupData, eval=TRUE,echo=TRUE, message=FALSE---------------------------------------------------------------------------------------------------------
library("gpuR")
library("clrng")
library("gpuBatchMatrix")
library('geostatsp')


## ----setupcoords, eval=TRUE--------------------------------------------------------------------------------------------------------------------------------
NlocalCache = 1000
Nglobal = c(128, 64, 2)
Nlocal = c(4, 2, 2)
theType = "double"


## ----swissRainBoundary, echo=TRUE--------------------------------------------------------------------------------------------------------------------------
data("swissRain", package="geostatsp")
myRaster = geostatsp::squareRaster(swissBorder, 80)
myRaster


## ----setup paramsBatch, eval=TRUE,echo=FALSE---------------------------------------------------------------------------------------------------------------
params = 
rbind(c(shape=1.25, range=50*1000, variance = 1.5, nugget = 0,anisoRatio = 1, anisoAngleRadians = 0), 
c(shape=2.15, range=60*1000, variance = 2, nugget = 0, anisoRatio = 4, anisoAngleRadians = pi/7),
c(shape=0.6, range=30*1000, variance = 2, nugget = 0, anisoRatio = 2, anisoAngleRadians = pi/7),
c(shape=3, range=30*1000, variance = 2, nugget = 0, anisoRatio = 2, anisoAngleRadians = pi/7)
)
#c(shape=2.15, range=40*1000, variance = 2, nugget = 0, anisoRatio = 2, anisoAngleRadians = pi/4))

## ----show parameterBatch,eval=TRUE-------------------------------------------------------------------------------------------------------------------------
params
#paramsGpu = gpuBatchMatrix::maternGpuParam(myParamsBatch, type=theType)


## ----r outputBatchCreate, eval=TRUE------------------------------------------------------------------------------------------------------------------------
#maternCov = vclMatrix(0, nrow(paramsGpu)*nrow(coordsGpu), 
#                         nrow(coordsGpu),type=theType)
maternCov = gpuBatchMatrix::maternBatch(
  params, myRaster,          
  Nglobal=c(128,64), Nlocal=c(16,4))
dim(maternCov)


## ----r cholBach, eval=TRUE---------------------------------------------------------------------------------------------------------------------------------
diagMat = gpuBatchMatrix::cholBatch(maternCov, 
                          Nglobal= c(128, 8), Nlocal= c(32, 8))


## ----r randomNormalsSim1, eval=TRUE, cache=FALSE-----------------------------------------------------------------------------------------------------------
streamsGpu <- createStreamsGpu(n=128*64)

## ----r randomNormalsSim2, eval=TRUE, cache=FALSE-----------------------------------------------------------------------------------------------------------
zmatGpu = clrng::rnorm(
  c(nrow(maternCov),2), streams=streamsGpu, 
                Nglobal=c(128,64),
                type = theType)


## ----r multLowerdiag_LDZ, eval=TRUE, cache=FALSE, results='hide'-------------------------------------------------------------------------------------------
simMat = vclMatrix(0, nrow(zmatGpu), ncol(zmatGpu), 
                   type = gpuR::typeof(zmatGpu))

gpuBatchMatrix::multiplyLowerDiagonalBatch(
  output=simMat, L=maternCov, 
  D=diagMat, B=zmatGpu,
  diagIsOne = TRUE,   
  transformD = "sqrt", 
  Nglobal, Nlocal, NlocalCache)


## ----r plotSimSetup, echo=TRUE, eval=TRUE, cache=FALSE-----------------------------------------------------------------------------------------------------
simRaster = raster::brick(myRaster, nl = ncol(simMat)*nrow(params))
values(simRaster) = as.vector(as.matrix(simMat))


## ----r theSubcapName, include=FALSE------------------------------------------------------------------------------------------------------------------------
names(simRaster) = apply(expand.grid('par',1:nrow(params), 
                                     'sim', 1:ncol(simMat)), 1, paste, collapse='')
theSubcap = gsub("par", "parameter ", names(simRaster))
theSubcap = gsub("sim", ", simuation ", theSubcap)


## ----maternplot, echo=FALSE, eval=TRUE, cache=FALSE, dev='png', out.width = '0.47\\textwidth', fig.height=3, fig.width=4.5, fig.cap="Simulated Gaussian random fields", fig.subcap = theSubcap, fig.align='default'----
myCol = mapmisc::colourScale(breaks=sort(unique(c(-6, -4, seq(-2, 2), 4, 6))), style='fixed', col='Spectral')
for(D in names(simRaster)) {
  mapmisc::map.new(simRaster)
  plot(simRaster[[D]], legend=FALSE, add=TRUE, col=myCol$col, breaks=myCol$breaks)
  plot(swissBorder, add=TRUE)
}
mapmisc::legendBreaks("right", myCol, inset=0)

