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


dyn.unload('../src/gmrfLik.so')
dyn.load('../src/gmrfLik.so')


obsCov = as.matrix(cbind(as.data.frame(theY), intercept=1, as.data.frame(thecov)))
Nrep = 2
obsCov = obsCov[,c(rep(1,Nrep), 2:ncol(obsCov))]

Snugget = c(Inf,10^seq(3,-3,len=81))


date()
old=geostatsp:::loglikGmrfOneRange(
    oneminusar=themodel['oneminusar'],
    Yvec=obsCov[,(1:2)], Xmat=obsCov[,-(1:2)], 
    NN=myNN, 
    propNugget=1/Snugget,
    shape=2,  boxcoxInterval=NULL,
    reml=TRUE,
    sumLogY = NULL,
    adjustEdges=FALSE,
    optimizer=FALSE)
date()


new= loglikGmrfOneRange(
    oneminusar=themodel['oneminusar'],
    Yvec=obsCov[,(1:2)], Xmat=obsCov[,-(1:2)], 
    NN=myNN, 
    propNugget=1/Snugget,
    shape=2,
    reml=TRUE,
    sumLogY = NULL,
    adjustEdges=FALSE,
    optimizer=FALSE)
date()
plot(log10(new['propNugget',1,-1]), new['m2logL.ml',1,-1])
lines(log10(old['propNugget',1,-1]), old['m2logL.ml',1,-1]-ncell(theY),col='blue')




dyn.unload('../src/gmrfLik.so')
dyn.load('../src/gmrfLik.so')
obsCov = as.matrix(cbind(as.data.frame(theY), intercept=1, as.data.frame(thecov)))
newBc = loglikGmrfOneRange(
    oneminusar=themodel['oneminusar'],
    Yvec=exp(as.data.frame(theY)[,1]),  
  Xmat=as.matrix(cbind( intercept=1, as.data.frame(thecov))), 
    NN=myNN, 
    propNugget=1/Snugget,
    fixBoxcox=FALSE,
    shape=2,
    reml=TRUE,
    sumLogY = NULL,
    adjustEdges=FALSE,
    optimizer=FALSE)


plot(log10(newBc['propNugget',1,-1]), newBc['m2logL.ml',1,-1])



 

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


