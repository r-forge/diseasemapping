library('geostatsp')



myraster = squareRaster(raster(extent(0,8000,0,6000), ncols=40,nrows=30))

myQ = maternGmrfPrec(myraster, 
    param=c(shape=2, oneminusar=0.1,
        conditionalVariance=100),
    adjustEdges=FALSE)

themodel = attributes(myQ)$param[['theo']]
maternShape = attributes(myQ)$param$theo['shape']

theU = geostatsp::RFsimulate(myraster, model=themodel, n=1)

thecov = myraster
values(thecov) = c(rep(0,ncell(thecov)/2),
    rep(4,ncell(thecov)/2))
names(thecov)='x'
beta.x=5
theY = theU + beta.x*thecov

obsCov = cbind(as.data.frame(theY), intercept=1, as.data.frame(thecov))

dyn.unload('../src/gmrfLik.so')
dyn.load('../src/gmrfLik.so')

Snugget = c(0,0.5,1)
stuff = .Call('gmrfLik',myQ, as.matrix(obsCov), Snugget)

Nxy = ncol(obsCov)
Nnugget = length(Snugget)

array(stuff[seq(1, Nnugget*Nxy^2)], dim=c(Nxy,Nxy,Nnugget),
  dimnames=list(colnames(obsCov), colnames(obsCov), Snugget))

matrix(stuff[-seq(1, Nnugget*Nxy^2)], ncol=4, 
    dimnames=list(Snugget, c('det','detreml','logL', 'logReL')))
 