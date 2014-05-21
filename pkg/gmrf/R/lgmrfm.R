lgmrfm = function(data,formula,covariates=NULL,
                  shape=1,nugget=0,
                  oneminusar=seq(0.01, 0.4,len=4), 
                  rangeInCells=NULL,
                  range=NULL,
                  NN=NNmat(data),mc.cores=1,
                  ...) {
  
  if(!is.null(covariates)){
    covariates = stackRasterList(covariates,template=data)
    data = stack(data,covariates)
  }
  dataOrig=data
  csiz = xres(dataOrig)
  data = as.data.frame(data)
  Yvec = data[,
              as.character(attributes(terms(formula))$variables)[2]
              ]
  Xmat = model.matrix(formula,data=data)
  
  thel = loglikGmrf(oneminusar=oneminusar,
                    Yvec=Yvec,Xmat=Xmat,
                    NN=NN,propNugget=nugget,
                    shape=shape,mc.cores=mc.cores,...)
  
  
  
  thesummary = summaryGmrfFit(thel, cellSize = csiz)
  thesummary$complete = list()
  thesummary$complete$array = thel
  thesummary$complete$oneminusar = thel['oneminusar',1,]
  thesummary$complete$propNugget = thel['propNugget',,1]
  thesummary$complete$range = csiz*sqrt(2*shape*(1-thesummary$complete$oneminusar)/thesummary$complete$oneminusar)
  thesummary$complete$rangeInCells = thesummary$complete$range / csiz
 
  
  mleparam = 
    thesummary$ml[ ,'mle']
  
  
  if(mleparam['propNugget']>0) {
    thesummary$predict = conditionalGmrf(
      param=mleparam,
      Yvec=Yvec,Xmat=Xmat,
      template=dataOrig, NN=NN,
      mc.cores=mc.cores,...)
  } else {
    thesummary$predict = raster::brick(
      raster(dataOrig), nl=2)
    names(thesummary$predict) = c('fixed','resid')
    values(thesummary$predict)=NA
    values(thesummary$predict[[1]]) =
      Xmat %*% mleparam[
        paste("beta.",colnames(Xmat),sep='')]
    values(thesummary$predict[['resid']]) =
      Yvec -values(thesummary$predict[['fixed']])
    
  }
  
  return(thesummary)
}

plotLgmrf = function(x,reml=FALSE){
  oldpar = par()
  par(mfrow=c(1,2))
  plotcols = c('mle','q0.025','q0.975')
  
  plot(x$profL$rangeInCells[[c('ml','reml')[1+reml] ]],type='o')
  abline(v=x[[c('ml','reml')[1+reml]]][
    'rangeInCells',plotcols
    ],col=c('red','orange','orange'))
  
  plot(x$profL$propNugget[,c(1,2+reml)],type='o')
  abline(v=x[[c('ml','reml')[1+reml]]][
    'propNugget',plotcols
    ],col=c('red','orange','orange'))
  par(oldpar)
}

plotLgmrf2 = function(x,reml=FALSE){
  oldpar = par()
  par(mfrow=c(1,2))
  plotcols = c('mle','q0.025','q0.975')
  
  plot(x$profL$rangeInCells[[c('ml','reml')[1+reml] ]],type='o')
  abline(v=x[[c('ml','reml')[1+reml]]][
    'rangeInCells',plotcols
    ],col=c('red','orange','orange'))
  
  plot(x$profL$propNugget[,c(1,2+reml)],type='o')
  abline(v=x[[c('ml','reml')[1+reml]]][
    'propNugget',plotcols
    ],col=c('red','orange','orange'))
  par(oldpar)
}