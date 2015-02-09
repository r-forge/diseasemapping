library('geostatsp')
data('swissRain')

Ncores = c(1,2)[1+(.Platform$OS.type=='unix')]



if(FALSE){
  dyn.unload("/home/patrick/workspace/diseasemapping/pkg/geostatsp/src/matern.so")
  dyn.load("/home/patrick/workspace/diseasemapping/pkg/geostatsp/src/matern.so")
}

sr2 = swissRain
sr2$elev = raster::extract(swissAltitude, sr2)
swissFit = likfitLgm(
    data=sr2, 
    formula=rain~ elev,
    param=c(range=10000,shape=1,nugget=0,boxcox=0.5,anisoRatio=5,anisoAngleDegrees=30),
    paramToEstimate = c("range",'anisoAngleDegrees','anisoRatio'),
    reml=TRUE
)

swissFit$param
swissFit$opt$boxcox
swissFit$opt$logL
swissFit$opt$total
swissFit$betaHat


sl = loglikLgm(
    swissFit$param[c('range','shape','boxcox', 'anisoRatio', 'anisoAngleRadians')],
    data=sr2, 
    formula=rain~ elev,
    reml=FALSE)  

attributes(sl)$totalVarHat

#sigSqHat =   (attributes(sl)$Lorig - attributes(sl)$totalVarHat) /100
sigSqHat = attributes(sl)$totalVarHat
sl1 = loglikLgm(
    c(attributes(sl)$param[c('boxcox','anisoRatio','anisoAngleRadians','shape', 'range')], variance=sigSqHat),
    data=sr2, 
    formula=rain~ elev,
    reml=FALSE)  

attributes(sl1)$totalVarHat

attributes(sl1)$Lorig
attributes(sl)$Lorig




2*attributes(sl1)$determinants[1]
2*attributes(sl)$determinants[1] + 100*log(sigSqHat)

2*attributes(sl1)$determinants[1] + attributes(sl1)$Lorig - attributes(sl1)$totalVarHat 

2*attributes(sl)$determinants[1] + 100*log(sigSqHat) + 100



as.numeric(sl) 
as.numeric(sl1)

attributes(sl)$Ltype
attributes(sl1)$Ltype

attributes(sl)$Lorig
attributes(sl1)$Lorig
attributes(sl)$totalVarHat
attributes(sl1)$totalVarHat
attributes(sl)$determinants
attributes(sl1)$determinants

attributes(sl)$param['variance']


mydet = c(
determinant(matern(sr2, param=attributes(sl)$param[c('anisoRatio','anisoAngleRadians','shape', 'range')]))$mod,
determinant(matern(sr2, param=
            c(attributes(sl)$param[c('anisoRatio','anisoAngleRadians','shape', 'range')],
                variance=sigSqHat)
    ))$mod
)
 c(mydet, mydet[1]+100*log(sigSqHat))

sigSqHat =   (attributes(sl)$Lorig - attributes(sl)$totalVarHat) /100
100*log(sigSqHat) + 2*attributes(sl)$determinants[1] + 100


2*attributes(sl1)$determinants[1]
2*attributes(sl)$determinants[1] + 100*log(sigSqHat)

attributes(sl)$betaHat
attributes(sl1)$betaHat

rbind(
    attributes(sl)$param,
attributes(sl1)$param)


as.numeric(sl)
attributes(sl)$twoLogJacobian
attributes(sl)$total
attributes(sl)$param
attributes(sl1)$varBetaHat
as.numeric(sl)+attributes(sl)$twoLogJacobian
 

x=profLlgm(swissFit, mc.cores=Ncores,
    boxcox=seq(0, 1 , len=6)
)
dots=list(    boxcox=seq(0, 1 , len=6));fit=swissFit


dots=list(		anisoAngleDegrees=seq(-20, 20 , len=6)+swissFit$param['anisoAngleDegrees'])


swissInf = informationLgm(swissFit)



if(!interactive()) 
  pdf("profLswissAngle.pdf")

plot(x[[1]],x[[2]], xlab=names(x)[1],
#		yaxt='n',
		ylab='log L',
		ylim=c(min(x[[2]]),x$maxLogL),
		type='n')
abline(h=x$breaks[-1],
		col=x$col,
		lwd=1.5)
axis(2,at=x$breaks,labels=x$prob,line=-1.2,tick=F,
		las=1,padj=1.2,hadj=0,col.axis='red')

abline(v=x$ciLong$par,
		lty=2,
		col=x$col[as.character(x$ciLong$prob)])


axis(1,at=x$ciLong$par,
		labels=x$ciLong$quantile,
		padj= -6,hadj=0.5, 
		tcl=0.5,cex.axis=0.8,
		col=NA,col.ticks='red',col.axis='red')

ciCols = grep("^ci", colnames(swissInf$summary),
		value=TRUE)

ciValues = unlist(swissInf$summary[intersect(names(x), rownames(swissInf$summary)),ciCols])
if(any(!is.na(ciValues)))
  axis(1,at=ciValues,
		labels=gsub("^ci","",ciCols),
		padj= 2,hadj=0.5, 
		tcl=-2,cex.axis=0.7,
		col=NA,col.ticks='blue',col.axis='blue')

lines(x[[1]],x[[2]])

if(!interactive()) 
  dev.off()


if(interactive()  | Sys.info()['user'] =='patrick') {
x2d=profLlgm(swissFit, mc.cores=Ncores,
		anisoAngleDegrees=seq(25, 35 , len=6),
		anisoRatio = exp(seq(log(4),log(16),len=8))
)
if(!interactive()) 
  pdf("profLswiss2d.pdf")
image(x2d[[1]],x2d[[2]],x2d[[3]],
		breaks=x2d$breaks,
		col=x2d$col,log='y',
		xlab=names(x2d)[1],
		ylab=names(x2d)[2])

thesevars = c("anisoAngleRadians","log(anisoRatio)")
thisV = swissInf$information[
		thesevars,thesevars]
thisMean= c(x2d$MLE["anisoAngleDegrees"],
		log(x2d$MLE['anisoRatio']))


if(requireNamespace("ellipse", quietly=TRUE)) {

for(D in x2d$prob[x2d$prob>0&x2d$prob<1]) {
	thisE = ellipse::ellipse(thisV, centre=thisMean,
			level=D)
  colnames(thisE) = names(thisMean)
	thisE = cbind(thisE,
			anisoRatioExp = exp(thisE[,"anisoRatio"]))
	lines(thisE[,"anisoAngleDegrees"],
			thisE[,"anisoRatioExp"],lwd=4)
	lines(thisE[,"anisoAngleDegrees"],
			thisE[,"anisoRatioExp"], col=x2d$col[as.character(D)],
			lwd=3)
}
}

points(x2d$MLE[1],x2d$MLE[2],pch=15) 

legend("topleft", fill=x2d$col, legend=x2d$prob[-length(x2d$prob)])

if(!interactive()) 
  dev.off()
}