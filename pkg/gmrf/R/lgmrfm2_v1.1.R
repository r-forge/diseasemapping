lgmrfm = function(data,formula,covariates=NULL,
                  shape=1,nugget=0,
                  oneminusar=seq(0.01, 0.4,len=4), 
                  rangeInCells=NULL,
                  range=NULL, reml = TRUE,
                  NN=NNmat(data),mc.cores=1,
                  ...) {
  
  if(!is.null(covariates)){
    covariates = stackRasterList(covariates,template=data)
    data = stack(data,covariates)
  }
  
  dataOrig=data
  csiz = xres(dataOrig)
  data = as.data.frame(data)
  Yvec = data[,as.character(attributes(terms(formula))$variables)[2]]
  Xmat = model.matrix(formula,data=data)
  thel = loglikGmrf(oneminusar=oneminusar,
                    Yvec=Yvec,Xmat=Xmat,
                    NN=NN,propNugget=nugget,
                    shape=shape,mc.cores=mc.cores,...)
  
  thesummary = list()
  if (reml){
  chooseLike = 'logL.reml'
  m2Like = 'm2logL.reml'
  }else{
  chooseLike = 'logL.ml'
  m2Like = 'm2logL.ml'
  }
  
  
  if (nugget == 0){
  thesummary$profL$propNugget = 0
  propNug = 0
  omar = thel['oneminusar',]
  rangeValue = csiz*sqrt(2*shape*(1-omar)/omar)
  rangeInCellsValue = rangeValue/csiz
  rangeLike = thel[chooseLike,]
  rangem2Like = thel[m2Like,]
  rangePack = cbind(omar, rangeValue, rangeInCellsValue, rangeLike, rangem2Like)
  rownames(rangePack) = NULL
  colnames(rangePack) = c('oneminusar', 'range', 'rangeInCells',chooseLike, m2Like)
  
  thesummary$profL$range = as.data.frame(rangePack)
  }
  else{
  # $propNugget
  propNug = thel['propNugget',,1]
  propLike =  thel[chooseLike,,1]
  propm2Like =  thel[m2Like,,1]
  propPack = cbind(propNug, propLike, propm2Like)
  rownames(propPack) = NULL
  colnames(propPack) = c('propNugget', chooseLike, m2Like)
  
  thesummary$profL$propNugget = as.data.frame(propPack)
  
  
  #$range
  omar = thel['oneminusar',1,]
  rangeValue = csiz*sqrt(2*shape*(1-omar)/omar)
  rangeInCellsValue = rangeValue/csiz
  rangeLike = thel[chooseLike,1,]
  rangem2Like = thel[m2Like,1,]
  rangePack = cbind(omar, rangeValue, rangeInCellsValue, rangeLike, rangem2Like)
  rownames(rangePack) = NULL
  colnames(rangePack) = c('oneminusar', 'range', 'rangeInCells',chooseLike, m2Like)
  
  thesummary$profL$range = as.data.frame(rangePack)
  }
  
  #$twoDim
  thesummary$profL$twoDim = list()
  thesummary$profL$twoDim$oneminusar = omar
  thesummary$profL$twoDim$range = rangeValue
  thesummary$profL$twoDim$rangeInCells = rangeInCellsValue
  thesummary$profL$twoDim$propNugget = propNug
  thesummary$profL$twoDim$array = thel
  
  


  thesummary$data = data
  thesummary$model$reml = reml
  thesummary$model$trend = formula
  if (reml){
  thesummary$summary = summaryGmrfFit(thel)$reml
  }else{
    thesummary$summary = summaryGmrfFit(thel)$ml  
  }

  thesummary$param = thesummary$summary[,'mle']
  return(thesummary)
}




# contourPlot = function(res, param = c('propNugget', 'range', 'rangeInCells')) { 
#   dseq = rev(c(0,0.5, 1,2,4,8,20))
#   
#   thecol = mapmisc::colourScale(res$profL$array['logL.ml',,],
#                                 breaks=max(res$profL$array['logL.ml',,]) - dseq,
#                                 col='RdYlGn',style='fixed',rev=TRUE)
#   count = grep(param, names(res$profL), value = F)[1]
#   para = res$profL[count]
#   plot(range(para),
#        range(res$profL$oneminusar),type='n',
#        xlab=param,ylab='1-ar')
#   .filled.contour(sort(para), 
#                   sort(res$profL$oneminusar),
#                   res$profL$array['logL.ml',,],
#                   col=thecol$col,levels=thecol$breaks)
# }

# 
# plotLgmrf2 = function(x,reml=FALSE){
#   oldpar = par()
#   par(mfrow=c(1,2))
#   plotcols = c('mle','q0.025','q0.975')
#   
#   plot(x$profL$rangeInCells[[c('ml','reml')[1+reml] ]],type='o')
#   abline(v=x[[c('ml','reml')[1+reml]]][
#     'rangeInCells',plotcols
#     ],col=c('red','orange','orange'))
#   
#   plot(x$profL$propNugget[,c(1,2+reml)],type='o')
#   abline(v=x[[c('ml','reml')[1+reml]]][
#     'propNugget',plotcols
#     ],col=c('red','orange','orange'))
#   par(oldpar)
# }

