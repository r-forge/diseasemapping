
# data is a raster.  grid is ignored
setMethod("lgm", 
    signature("formula", "Raster", "ANY", "ANY"),
    function(
        formula, 
        data,  
        grid=NULL,
        covariates=NULL, ...) {
      
      dataCov = gm.dataRaster(
          formula, data,
          grid=raster(data),
          covariates=covariates,
          buffer=0)
      
      callGeneric(formula, 
          data=dataCov$data, 
          grid=dataCov$grid, 
          covariates=dataCov$covariates, ...)
    }
)


setMethod("lgm", 
    signature("formula", "data.frame", "Raster", "data.frame"), 
    function(formula,
        data, grid,
				covariates=NULL,
        shape=1,boxcox=1,nugget=0,
				  expPred=FALSE, nuggetInPrediction=TRUE,
          reml = TRUE, mc.cores=1,
          fixBoxcox=TRUE,
          fixNugget = FALSE,
          ...)
      {
  NN=NNmat(grid)
  csiz = xres(grid)
        
	Yvar = all.vars(formula)[1]
	if(!length(grep('[[:digit:]]$', Yvar))) {
	  allYvar = grep(paste("^", Yvar, "[[:digit:]]*$",sep=""), names(data), value=TRUE)
  } else {
		allYvar = Yvar
	}
  Yvec = as.matrix(data[,allYvar, drop=FALSE])
	
	if(!Yvar %in% colnames(data)) {
		data[,Yvar] = 0
	}
  Xmat = model.matrix(formula, data)		  

	if(nrow(Xmat) != ncell(grid))
		warning("dimensions of data and grid are not compatible")
  
  if(!fixNugget & (length(nugget)<2))
    nugget = NULL
  
  thel = loglikGmrf(Yvec=Yvec,Xmat=Xmat,
                    NN=NN, 
                    propNugget=nugget,
                    boxcox=boxcox, fixBoxcox=fixBoxcox,
                    shape=shape,mc.cores=mc.cores,
                    reml=reml, ...)
  mle = thel$mle   
  lArray = thel$mlArray
  lMat = thel$mlMat
  

  if (reml){
    chooseLike = 'logLreml'
		varInMle = c('tausqReml', 'varReml')
  }else{
    chooseLike = 'logLml'
		varInMle = c('tausqMl', 'varMl')
  }
  

  res = list(param = drop(mle))
  if(!is.null(lArray)){
    res$array = lArray
    res$profL = list()
    # nugget
    if(dim(lArray)[2]>1){ # have nugget
      best = apply(lArray[,,chooseLike,,drop=FALSE],
          2, which.max) 
      best=arrayInd(best, dim(lArray)[-(2:3)])
      res$profL$propNugget = NULL
      for(D in 1:nrow(best)){
        res$profL$propNugget = rbind(
            res$profL$propNugget,
            lArray[best[1], D, 
                c('propNugget',chooseLike),
                best[2]])
      }
    } # end have nugget
     
    if(dim(lArray)[4]>1){ # have oneminusar
      best = apply(lArray[,,chooseLike,,drop=FALSE],
          4, which.max) 
      best=arrayInd(best, dim(lArray)[-(3:4)])
      
      res$profL$oneminusar = NULL
      for(D in 1:nrow(best)){
        res$profL$oneminusar = rbind(
            res$profL$oneminusar,
            lArray[best[1], best[2], 
                c('oneminusar',chooseLike,'range'),
                D])
      }
    } # end have oneminusar
    res$profL$range = res$profL$oneminusar[,c(3,2)]
    res$profL$oneminusar = res$profL$oneminusar[,c(1,2)]
    
    if(all(dim(lArray)[c(2,4)]>1)){ # have both
      x=lArray[1,,'propNugget',1]
      orderx = order(x)
      y = lArray[1,1,'range',]
      ordery=order(y)
      
      res$profL$twoDim = list(
          x=lArray[1,orderx,'propNugget',1],
          y=lArray[1,1,'range',ordery],
          z=apply(
          lArray[,,chooseLike,,drop=FALSE],
          c(2,4), max, na.rm=TRUE)[orderx,ordery],
      oneminusar=lArray[1,1,'oneminusar',ordery]      
    )
 
    } # end have both
  } # end lArray not null

  
  res$data = data
  res$model$reml = reml
  res$model$trend = formula
 
  # summary table
	covInMle = grep("Se$", rownames(mle), value=TRUE)

  scovariates = gsub(
      'Se$','', covInMle
  )
  
  srownames = c('sdNugget','sdSpatial','range','shape')

  scolnames = c("estimate", "stdErr", "ci0.005", "ci0.995", "ci0.025", "ci0.975", 
      "ci0.05", "ci0.95", "ci0.1", "ci0.9", "pval", "Estimated")
	ress = list()
	for(D in 1:ncol(mle)){
    ress[[D]] = as.data.frame(
      matrix(
          NA,
          length(scovariates) + length(srownames),
          length(scolnames),
          dimnames = list(
              c(scovariates, srownames),
              scolnames
              )
          )
      )
   ress[[D]][c('sdNugget','sdSpatial','range','shape'),'estimate'] = 
       c(sqrt(mle[c('tausq','sigmasq'),D]),
        mle[c('range','optimalShape'),D]   
       )
   ress[[D]][c('sdNugget','sdSpatial','range','shape'),'Estimated'] =
       c(fixNugget, TRUE, TRUE, FALSE)
   ress[[D]][scovariates,'Estimated']  = TRUE   
   ress[[D]][scovariates,'estimate']  = mle[covInMle,D]   
   ress[[D]][scovariates,'stdErr']  = mle[covInMle,D]   
 } # for D
       
  if(length(ress)==1) ress = ress[[1]]
	  res$summary = ress

	
  return(res)
  

}
)

