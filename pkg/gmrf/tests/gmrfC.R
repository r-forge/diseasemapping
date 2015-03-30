library('geostatsp')



myraster = squareRaster(raster(extent(0,8000,0,6000), ncols=20,nrows=30))
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
theY = theU + beta.x*thecov

obsCov = as.matrix(cbind(as.data.frame(theY), intercept=1, as.data.frame(thecov)))
Nrep = 2
obsCov = obsCov[,c(rep(1,Nrep), 2:ncol(obsCov))]

YrepAdd = rep(0, Nrep)

Snugget = c(Inf,10^seq(3,-3))
Nnugget= length(Snugget)
Nxy = ncol(obsCov)

obsCov = as.matrix(obsCov)
dyn.unload('../src/gmrfLik.so')
dyn.load('../src/gmrfLik.so')

stuff1 = .Call('gmrfLik',myQ, obsCov[,-2], Snugget, as.double(YrepAdd[-2]))

ssq1=array(stuff1[seq(1, Nnugget*(Nxy-1)^2)], 
    dim=c(Nxy-1,Nxy-1,Nnugget),
    dimnames=list(colnames(obsCov)[-2], colnames(obsCov)[-2], Snugget))


stuff2 = .Call('gmrfLik',myQ, obsCov, Snugget, as.double(YrepAdd))

ssq2=array(stuff2[seq(1, Nnugget*Nxy^2)], 
    dim=c(Nxy,Nxy,Nnugget),
    dimnames=list(colnames(obsCov), colnames(obsCov), Snugget))



ssq1[,,1]
ssq2[,,1]

ml = array(
    stuff2[-seq(1, Nnugget*Nxy^2)], 
    dim=c(Nrep, Nnugget, 6),
    dimnames=list(names(YrepAdd), Snugget, 
        c('det','detreml','logL', 'logReL', 'varMl', 'varReml')))

ml = cbind(ssq=ssq[1,1,], detSum = - ml[1,'det']+ml[,'det'], ml)





xisqTausq = 10
V =  xisqTausq*solve(myQ) + diag(ncol(myQ))
Vchol = t(chol(V))
C = myQ + xisqTausq*diag(ncol(myQ))
Cchol = t(chol(C))
2*determinant(Vchol, modulus=TRUE)$mod
2*determinant(Cchol, modulus=TRUE)$mod -2*determinant(Qchol, modulus=TRUE)$mod
2*determinant(Cchol, modulus=TRUE)$mod 
2*determinant(Qchol, modulus=TRUE)$mod

# hard way
yvy = crossprod(solve(Vchol, obsCov))
#yvy = t(obsCov) %*% myQ %*% obsCov
betahat= solve(yvy[-1,-1])%*% yvy[-1,1]
ssqHard = as.vector(yvy[1,1] - yvy[-1,1]%*%betahat)
ssqHard

# easy way
xy = crossprod(obsCov)
Lxy = solve(Cchol, obsCov)
yvyE = xy - xisqTausq*crossprod(Lxy)
betahatE = as.matrix(solve(yvyE[-1,-1])%*%yvyE[-1,1])

as.matrix(yvy)
as.matrix(yvyE)
ssq[,,as.character(xisqTausq)]
betahatE
betahat

theldl = Cholesky(C, LDL=FALSE, perm=TRUE)
obsCovRot = solve(theldl,as.matrix(obsCov), system='P')
crossprod(solve(theldl, obsCovRot,
        system='L'))
xisqTausq*crossprod(Lxy)

xy = crossprod(obsCov)
Lxy = solve(Cchol, xy)
yvyE = xy - xisqTausq*crossprod(Lxy)

solve(yvy[-1,-1])

Nxy = ncol(obsCov)
Nnugget = length(Snugget)


betahat
yvy[-1,1]
ssq[1,-1,1]
ssq[-1,1,1]
sum(ssq[1,-1,1]*ssq[-1,1,1])

as.matrix(yvy)
ssq[,,1]
solve(yvy[-1,-1])

Qchol = Cholesky(myQ,LDL=FALSE)

2*determinant(Qchol, modulus=TRUE)$mod



fromgp = loglikGmrfGivenQ(
   1/xisqTausq,
    obsCov[,1,drop=FALSE],
    obsCov[,-1],
    myQ, 
    Qchol=NULL, 
    detQ=NULL,
    boxcoxInterval=NULL,
    reml=FALSE
)

fromgp[c(1:2),]
fromgp[13:14,]
fromgp[5:6,]*(nrow(obsCov)-c(0,2))
ssq[,,as.character(xisqTausq)]
ssq[2,3,] = ssq[3,2,]
solve(ssq[-1,-1,as.character(xisqTausq)])



res=NULL
for(Dnugget in Snugget)
  res = cbind(res, 
      geostatsp:::loglikGmrfGivenQ(
    1/Dnugget,
    as.matrix(obsCov[,1]),as.matrix(obsCov[,-1]),
    myQ, 
    Qchol=Qchol, 
    detQ=NULL,
    boxcoxInterval=NULL,
    reml=FALSE
  )[c(1,2,13,14),])

rbind(Snugget,res)
ssq[-1,1,]




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