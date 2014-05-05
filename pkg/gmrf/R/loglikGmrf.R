objfunM = function(oparam,distVec,sqrtVar,
		range=NULL,shape=NULL){
	
	if(!is.null(range))
		oparam['range'] = as.vector(range)
	if(!is.null(shape))
		oparam['shape'] = as.vector(shape)
	if(!any(names(oparam)=='variance')){	
		oparam['variance']=1
	}
	if(!any(names(oparam)=='nugget')){	
		oparam['nugget']=1
	}
	
	theM =
			geostatsp::matern(
			x=distVec, 
			param=oparam[c('range','shape','variance')]
	)
thezeros = which(distVec<=0)
	theM[thezeros] = theM[thezeros] + oparam['nugget']
	
	sum( (sqrt(theM) - sqrtVar)^2)
}

loglikGmrfGivenQ = function(
		propNugget, Yp, YL, 
		Xp, XL,Qchol, detQ, 
		QLL=NULL,cholQLL=NULL,
		empirical=NULL) {
	

	# propNugget =    tau^2/xi^2
	
	
	if(propNugget>0){


		# Q = L L', QLL = L'L
		# cholIcQ = chol of L'L + I xisq/tausq
		if(is.null(QLL)) {
			QLL = crossprod(expand(Qchol)$L)

			cholIcQ = Cholesky(QLL, LDL=FALSE,
					Imult=1/propNugget)
			YL = solve(cholIcQ, YL,
					system='P')	
			XL = solve(cholIcQ, XL,
					system='P')	
		} else {
			cholIcQ = update(cholQLL, QLL,mult=1/propnugget)
		}

		
		Ybreve = solve(cholIcQ, 
				YL,	system='L')
		Xbreve = solve(cholIcQ, 
				 XL, system='L')	
		# should divide the above by sqrt(nuggetSigsq)
		
		XprecXinv = solve(Matrix::crossprod(Xbreve,Xbreve))
		
		betaHat = as.vector(
				XprecXinv %*% 
						Matrix::crossprod(Xbreve,Ybreve) 
		)
		names(betaHat) = colnames(XL)
		
		XprecXinv=XprecXinv*nuggetSigsq
		
		residsL =  as.vector(Ybreve -  Xbreve %*% betaHat )/ 
				sqrt(nuggetSigsq)
		
		
		logDetVar = as.numeric(
				length(YL) * log(nuggetSigsq) + 
				2*determinant(cholIcQ,logarithm=TRUE)$modulus -
						detQ 
		)
		

		
	} else { # no nugget 
		nuggetSigsq = 0

		
		XprecX = Matrix::crossprod(XL,XL)
		XprecXinv = solve(XprecX)
		
		betaHat = as.numeric(XprecXinv %*% 
						Matrix::crossprod(XL , YL))
		names(betaHat) = colnames(XL)

		
		# residsL ~ N(0,sigsq*I)
		residsL  = as.vector(YL - XL %*% betaHat)
 		
		logDetVar=-detQ
		# matrix operations
		# Qchol, cholCovMat = Matrix::chol(covMat)
		# XL, cholCovInvX = Matrix::solve(cholCovMat, covariates)
		# XprecX, cholCovInvXcross = Matrix::crossprod(cholCovInvX)
		# XprecXinv, cholCovInvXcrossInv = Matrix::solve(cholCovInvXcross)
	}

 
	
	rpr = as.numeric(crossprod(residsL,residsL))

	N=c(ml=0,reml=-ncol(Xp))+length(Yp)
	# estimate of the constant 
	constHat = rpr/N


	m2logL = N*log(constHat) + logDetVar + 
			c(ml=0,reml=determinant(
							XprecXinv,logarithm=TRUE)$modulus
	)
	names(m2logL) = paste("m2logL.",names(m2logL),sep='')


	variances = rep(c(
					sigmasq = 1,
					tausq = nuggetSigsq),2)*
		rep(constHat,c(2,2)) 
	names(variances) = paste(names(variances),
			rep(names(N),c(2,2)),sep='.')

	logL = -m2logL/2
	names(logL) = gsub("^m2","",names(logL))




	
	sebeta = diag(XprecXinv)
	sebeta = rep(sebeta,2)*
			rep(constHat,
				rep(length(sebeta),2))
	
#	print(sebeta)
	
	names(sebeta)= paste('se',rep(names(betaHat),2), 
			names(sebeta),sep='.')
	
	varbetahat = abind::abind(
			as.matrix(XprecXinv)/constHat[1],
			as.matrix(XprecXinv)/constHat[2],
			along=3)
			
	dimnames(varbetahat) = 
					list(names(betaHat),names(betaHat),
							names(constHat))

	if(!is.null(empirical)){

		empirical$yMult = empirical$yMult * variances['sigmasq.ml']	
		empirical[empirical$x<=0,'yMult'] = 
				empirical[empirical$x<=0,'yMult'] + 
				variances['tausq.ml']
		startparam = c(
				shape=2,
				range=2,
				variance=as.vector(variances['sigmasq.ml'])
		)
		if(variances['tausq.ml']>0)
				startparam = c(startparam,
						nugget=as.vector(variances['tausq.ml'])				
				)

			postfit=optim(fn=objfunM,par=startparam,
				lower=c(shape=0.1,range=0.1,
						variance=0,nugget=0)[names(startparam)],
				upper=c(shape=3,range=12,
						variance=Inf,nugget=Inf)[names(startparam)],
				distVec=empirical$x, 
				sqrtVar = sqrt(empirical$yMult))$par
		
		
		
		names(postfit) = paste(names(postfit), ".postfit",sep='')
		
	} else {
		postfit=NULL
	}		
			
	result =c(m2logL, logL,
			variances,
			propNugget=propNugget,			
			beta=betaHat, 
			sebeta,
			logDetVar=as.numeric(logDetVar),
			rpr=rpr,postfit)
	
	
	attributes(result)$varbetahat = varbetahat
	
	result

}

loglikGmrfOneRange = function(
		oneminusar=NULL, rangeInCells=NULL,
		Yvec, Xmat, NN, propNugget=0,
		shape=1,
		adjustEdges=FALSE,adjustParam=FALSE,
		adjustShape=FALSE,adjustMarginalVariance=FALSE) {
	

	
	if(is.null(oneminusar)){

		NN = maternGmrfPrec(NN,
				param=c(shape=shape,
						range=as.vector(rangeInCells[1]),
						conditionalVariance=1),
				adjustEdges=adjustEdges,adjustParam=adjustParam,
				adjustShape=adjustShape,
				adjustMarginalVariance=adjustMarginalVariance)
		oneminusar=NA
	} else {
 		NN =  maternGmrfPrec(NN,
				param=c(shape=as.vector(shape),
						oneminusar=as.vector(oneminusar[1]),
						conditionalVariance=1),
				adjustEdges=adjustEdges,adjustParam=adjustParam,
				adjustShape=adjustShape,
				adjustMarginalVariance=adjustMarginalVariance)

		rangeInCells = 	
			attributes(NN)$param$target['rangeInCells']
		shape=attributes(NN)$param$target['shape']
		
				
	}		

	cholPrec = Cholesky(NN,LDL=FALSE)
	LofQ = expand(cholPrec)$L
	
	theDet = 2*as.numeric(determinant(cholPrec,
					logarithm=TRUE)$modulus)
	
	Xp = solve(cholPrec,Xmat,system="P")
	Yp = solve(cholPrec,Yvec,system="P")


	
	XL = Matrix::crossprod(LofQ,Xp)
	YL = as.vector(
			Matrix::crossprod(LofQ,Yp)
	)

	argList = list(Yp=Yp,YL=YL,
			Xp=Xp,XL=XL,
			Qchol=cholPrec,
			detQ=theDet,
			empirical=attributes(NN)$param$empirical)

	if(!all(propNugget==0)) {
		# calculate t(LofQ) LofQ
		pRev = as(ncol(LofQ):1, "pMatrix")		
		lQLL =  as( t(LofQ %*% pRev) %*% pRev,'dtCMatrix')
		QLL = tcrossprod(lQLL)
		QLL = forceSymmetric(QLL)

		
		cholQLL = Cholesky(QLL,LDL=F,perm=T)
		
		# QLL = P' L  L' P 
		# QLLorig = Prev P' L   L' P Prev'
	
		ptwice =   as(expand(cholQLL)$P,'sparseMatrix') %*% 
				t(as(pRev,'sparseMatrix'))
		
		ptwice2 = as(ptwice, 'pMatrix')
		cholQLL@perm = as.integer(ptwice2@perm-1)

		argList$cholQLL = cholQLL
		if(FALSE){
			#prove that cholQLL is cholesky of QLLorig
			QLLorig = crossprod(LofQ)
			eye = solve(cholQLL,QLLorig)
			range(diag(eye))
			range(eye[which(lower.tri(eye,diag=FALSE))])
		}

		
		argList$YL = solve(argList$cholQLL, argList$YL,
				system='P')	
		argList$XL = solve(argList$cholQLL, argList$XL,
				system='P')	
	}
	
# variance of Q inverse
theraster = attributes(NN)$raster
Nx=ncol(theraster);Ny=nrow(theraster)

midcellCoord = c(round(Nx*.5),round(Ny*0.5)) # the middle cell
midcell = c(Nx*(Ny-midcellCoord[2]) + midcellCoord[1])
midVec = sparseMatrix(midcell,1,x=1,
		dims=c(ncol(NN),1))

Qcentre = solve(cholPrec,
		midVec)[midcell]


	if(length(propNugget)==1) {
		
		argList$propNugget=propNugget

		res = do.call(loglikGmrfGivenQ,argList)

		res = c(res,
				rangeInCells=
						as.numeric(rangeInCells),
				oneminusar=as.numeric(oneminusar),
				shape=as.numeric(shape),
				Qcentre=as.numeric(Qcentre)
			)
		
	} else {

		res = mapply(loglikGmrfGivenQ,
			propNugget=propNugget,			
			MoreArgs=argList,SIMPLIFY=TRUE
		)
		res = rbind(res,
			rangeInCells=as.numeric(rangeInCells),
			oneminusar=as.numeric(oneminusar),
			shape=as.numeric(shape),
			Qcentre=as.numeric(Qcentre)
		)
	
}

attributes(res)$Qinfo = attributes(NN)$param

res

}



loglikGmrf = function(
		oneminusar=NULL, rangeInCells=NULL,
		Yvec, Xmat, NN, propNugget=0,
		shape=1,
		adjustEdges=FALSE,adjustParam=FALSE,
		adjustShape=FALSE, 
		adjustMarginalVariance=FALSE,
		mc.cores=1) {

	
		argList = list(Yvec=Yvec, 
			Xmat=Xmat, NN=NN,
			propNugget=propNugget,
			shape=shape,
			adjustEdges=adjustEdges,
			adjustParam=adjustParam,
			adjustShape=adjustShape,
			adjustMarginalVariance=adjustMarginalVariance
		)

 		if(mc.cores>1) {
			myapply= function(...){
				parallel::mcmapply(...,
						mc.cores=mc.cores)
			}
		} else {
			myapply= mapply
				}
		
		if(!is.null(oneminusar)){
			res=myapply(loglikGmrfOneRange,
				oneminusar=oneminusar,	
				MoreArgs=argList,
				SIMPLIFY=FALSE
				)
 				names(res) = 
						paste("oneminusar=",oneminusar,sep="")
		} else {
			res=myapply(loglikGmrfOneRange,
					rangeInCells=rangeInCells,
					MoreArgs=argList,SIMPLIFY=FALSE
			)
			names(res) = 
					paste("rangeInCells=",rangeInCells,sep="")
			
		}
		if(class(res[[1]])!= 'character'){
			Qinfo = attributes(res[[1]])$Qinfo

			res= simplify2array(res)
			
			attributes(res)$Qinfo = Qinfo
			
		}
		res
}

 


summaryGmrfFit= function(x) {
	UseMethod("summaryGmrfFit")	
}
 
summaryGmrfFit.array = function(x) {
	

	x2 = aperm(x,c(3,2,1))
	x2 = matrix(c(x2), ncol=dim(x2)[3])
	colnames(x2) = dimnames(x)[[1]]
	res = summaryGmrfFit.matrix(x2,npar=2)
	
	res$profL = list()
	res$profL$propNugget = cbind(propNugget=x['propNugget',,1],
			logL.ml=apply(x['logL.ml',,],1,max),
			logL.reml = apply(x['logL.reml',,],1,max)
			)
whichL = apply(x['logL.ml',,],2,which.max)

			res$profL$oneminusar = cbind(
			oneminusar = x['oneminusar',1,],
	logL.ml=apply(x['logL.ml',,],2,max),
			logL.reml=apply(x['logL.reml',,],2,max)
	)
	
	if(all(
			c('range.postfit','shape.postfit') %in%
					dimnames(x)[[1]])) {
		res$profL$oneminusar  = cbind(
				res$profL$oneminusar,
	range.postit = diag(
			x['range.postfit',whichL
					,]),
	shape.postit = diag(
			x['shape.postfit',whichL
					,])
	)
	}
	
	
	res$profL$rangeInCells = list(
			ml = cbind(
				rangeInCells=	x['rangeInCells',1, ],
				logL=	apply(x['logL.ml',,],2,max)
			),
			reml = cbind(
					rangeInCells=x['rangeInCells',1, ],
					logL=apply(x['logL.reml',,],2,max)
			)
	)
			
	res	
	
}

summaryGmrfFit.matrix = function(x,npar=1) {
	
	if(any(rownames(x)=='logL.ml')){
		x = t(x)
	}
result = list()
someL = c('ml','reml')
	for(D1 in someL) {
	D=paste('logL.',D1,sep='')
	MLErow = x[which.max(x[ ,D]),]	
	thebetas = grep("^beta",names(MLErow), value=TRUE)
	betaMLE = MLErow[thebetas]	

	thebetas2 = gsub("^beta\\.","",thebetas)
	betaSE = MLErow[paste("se.",thebetas2,".ml",sep="")]
	
	betamat = cbind(mle=betaMLE,
			se=betaSE,
			q0.025=betaMLE - 1.96*betaSE,
			q0.975=betaMLE + 1.96*betaSE,
			pval = 2*pnorm(
					abs(betaMLE)/
							betaSE,
					lower=FALSE)	
	)
	
	withinCI = which(-2*(x[,D] - 
					MLErow[D]) < qchisq(0.95, 2))
	
xsub = x[,grep("^se\\.",colnames(x),invert=TRUE),drop=FALSE]
	parCI= t(apply(
					xsub[withinCI, ,drop=FALSE],
					2,range))
	colnames(parCI) = c('q0.025', 'q0.975')
	
	parMat = cbind(mle=MLErow[colnames(xsub)],parCI,se=NA,pval=NA)
	parMat[rownames(betamat),colnames(betamat)] = betamat
	
	notD = someL[someL != D1]
	notD = grep(
			paste("\\.",notD, "$",sep=''),
			rownames(parMat), 
			invert=TRUE,value=TRUE)
	
	parMat = parMat[notD,]
	rownames(parMat) = gsub(paste("\\.", D1, "$",sep=""),
			"",rownames(parMat))
	
	result[[D1]] = parMat
	}
	return(result)		
}
 
