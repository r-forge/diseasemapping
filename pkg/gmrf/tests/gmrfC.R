library('geostatsp')



myraster = squareRaster(raster(extent(0,8000,0,6000), ncols=80,nrows=120))
myNN = NNmat(myraster)
myQ = maternGmrfPrec(myNN, 
    param=c(shape=2, oneminusar=0.1,
        conditionalVariance=1),
    adjustEdges=FALSE)

themodel = attributes(myQ)$param[['theo']]
maternShape = attributes(myQ)$param$theo['shape']

set.seed(0)
theU = RFsimulate(x=myraster, model=themodel, n=1)

thecov = myraster
values(thecov) = c(rep(0,ncell(thecov)/2),
    rep(4,ncell(thecov)/2))
names(thecov)='x'
beta.x=5
theY = theU + beta.x*thecov+4
values(theY) = rnorm(ncell(theY), mean=values(theY), sd=sqrt(0.1))


dyn.unload('../pkg/geostatspDevel/src/gmrfLik.so')
dyn.load('../pkg/geostatspDevel/src/gmrfLik.so')


Yvec=exp(as.data.frame(theY)[,1])
Xmat=as.matrix(cbind( intercept=1, as.data.frame(thecov)))
NN=myNN
propNugget=NULL
#propNugget=c(0, 0.25, 0.5)
fixBoxcox=FALSE
shape=2
reml=FALSE
adjustEdges=FALSE
oneminusar = NULL
#oneminusar = c(0.1, 0.2, 0.3)
boxcox=1
fixBoxcox=FALSE
seqBoxcox=outer(c(0.01, 0.02), c(-1,1))
jacobian = 0
control=list(
    oneminusar = c(lower=1e-3, upper=0.8, tol= .Machine$double.eps^0.25),
    propNugget = c(lower=1e-3, upper=1e3, tol= .Machine$double.eps^0.25)
)


obsCov = as.matrix(cbind(as.data.frame(theY), intercept=1, as.data.frame(thecov)))
Nrep = 2
obsCov = obsCov[,c(rep(1,Nrep), 2:ncol(obsCov))]

Snugget = c(Inf,10^seq(3,-3,len=201))


date()
old=geostatsp:::loglikGmrfOneRange(
    oneminusar=themodel['oneminusar'],
    Yvec=obsCov[,(1:2)], Xmat=obsCov[,-(1:2)], 
    NN=myNN, 
    propNugget=1/Snugget,
    shape=2,  boxcoxInterval=NULL,
    reml=TRUE,
    sumLogY = NULL,
    adjustEdges=FALSE)
date()


new= loglikGmrfOneRange(
    oneminusar=themodel['oneminusar'],
    Yvec=obsCov[,(1:2)], Xmat=obsCov[,-(1:2)], 
    NN=myNN, 
    propNugget=1/Snugget,
    shape=2,
    reml=TRUE,
    adjustEdges=FALSE)
date()
plot(log10(new['propNugget',1,-1]), new['m2logL.ml',1,-1])
lines(log10(old['propNugget',1,-1]), old['m2logL.ml',1,-1]-ncell(theY),col='blue')




dyn.unload('../src/gmrfLik.so')
dyn.load('../src/gmrfLik.so')
 newBc = loglikGmrfOneRange(
    oneminusar=themodel['oneminusar'],
    Yvec=exp(as.data.frame(theY)[,1]),  
  Xmat=as.matrix(cbind( intercept=1, as.data.frame(thecov))), 
    NN=myNN, 
    propNugget=1/Snugget,
    fixBoxcox=FALSE,
    shape=2,
    reml=FALSE,
    adjustEdges=FALSE)


plot((newBc[,'propNugget']), newBc[,'m2logL.ml'], log='x')
abline(v=newBc[1,'propNugget'], col='red')



dyn.unload('../src/gmrfLik.so')
dyn.load('../src/gmrfLik.so')

newOpt = loglikGmrfOneRange(
    oneminusar=themodel['oneminusar'],
    Yvec=exp(as.data.frame(theY)[,1]),  
    Xmat=as.matrix(cbind( intercept=1, as.data.frame(thecov))), 
    NN=myNN, 
    propNugget=NULL,
    fixBoxcox=FALSE,
    shape=2,
    reml=FALSE,
    adjustEdges=FALSE)

xcol = 'xisqTausq'
nearMin = which(newBc[,'m2logL.ml'] - min(newBc[,'m2logL.ml'] ) < 12)
plot((newBc[nearMin,xcol]), newBc[nearMin,'m2logL.ml'], 
    log='x', col='red', 
    ylim=range(
        c(newBc[nearMin,'m2logL.ml'],min(newOpt[,'m2logL.ml']))
  ))
abline(v=newBc[1,xcol], col='red')
points(newOpt[,xcol], newOpt[,'m2logL.ml'],
    col='blue')
abline(v=newOpt[1,xcol], col='blue', lty=3)

V =  xisqTausq[2]*solve(myQ) + diag(ncol(myQ))
Vchol = t(chol(V))


C = myQ + xisqTausq[2]*diag(ncol(myQ))
Cchol = t(chol(C))
Qchol = Cholesky(Q,LDL=FALSE, perm=FALSE)
2*determinant(Vchol, modulus=TRUE)$mod
2*determinant(Cchol, modulus=TRUE)$mod -2*determinant(Qchol, modulus=TRUE)$mod
2*determinant(Cchol, modulus=TRUE)$mod 
2*determinant(Qchol, modulus=TRUE)$mod


yvy = crossprod(solve(Vchol, obsCov))
betahat=  solve(yvy[-(1:Ny),-(1:Ny)])%*% yvy[-(1:Ny),(1:Ny)]

yvytop = yvy[1:Ny,1:Ny] - yvy[(1:Ny),-(1:Ny)]%*% betahat

Nobs * log(diag(yvytop)) - Nobs * log(Nobs) + 2*determinant(Vchol, modulus=TRUE)$mod +YrepAdd[2]


ml[,2,]

t(obsCov)  %*% Q %*%  obsCov
yvy


# hard way
yvy = crossprod(solve(Cchol, obsCov))
yvy
#yvy = t(obsCov) %*% myQ %*% obsCov

ssqHard =    yvy[(1:Ny), (1:Ny)] - t(yvy) %*%t(betahat) 
ssqHard 
ssq[,,2]


