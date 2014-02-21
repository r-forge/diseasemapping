loglikGmrfGivenQ = function(
		propNugget, Yp, YL, 
		Xp, XL,Qchol, detQ, 
		QLL=NULL,cholQLL=NULL) {
	
	if(propNugget>0){

		# propNugget = tausq / (tausq + sigsq)
		nuggetSigsq = 1/( (1/propNugget) -1)
		# nuggetSigsq = tausq/sigsq
		
		
		if(is.null(QLL)) {
			LofQ = expand(Qchol)$L
			QLL = forceSymmetric(crossprod(LofQ,LofQ))
			diag(QLL) = diag(QLL) + 1/nuggetSigsq
			cholIcQ = Cholesky(QLL, LDL=FALSE)
			YL = solve(cholIcQ, YL,
					system='P')	
			XL = solve(cholIcQ, XL,
					system='P')	
		} else {
			# IcQ = (Q + I/c) c 
			cholIcQ = update(cholQLL, QLL,mult=1/nuggetSigsq)
		}
		# should multiply cholIcQ by sqrt(nuggetSigsq)

		
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
		names(betaHat) = colnames(Xmat)
		
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
		names(betaHat) = colnames(Xmat)
		
		# resids ~ N(0,sigsq*prec^(-1))
		resids = Yp - as.vector(Xp %*% betaHat)
		residsL = as.numeric(Matrix::crossprod(
				resids, 
				expand(Qchol)$L))
		
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
		rep(rpr/N,c(2,2)) 
	names(variances) = paste(names(variances),
			rep(names(N),c(2,2)),sep='.')

	logL = -m2logL/2
	names(logL) = gsub("^m2","",names(logL))

	
	varbetahat = as.matrix(XprecXinv)
	sebeta = diag(varbetahat)
	sebeta = rep(sebeta,2)*rep(constHat,rep(length(sebeta),2))
	names(sebeta)= paste('se',colnames(Xmat), names(sebeta),sep='.')
	
	varbetahat = abind::abind(
			varbetahat/constHat[1],
			varbetahat/constHat[2],
			along=3)
			
	dimnames(varbetahat) = 
					list(names(betaHat),names(betaHat),
							names(constHat))

			
	result =c(m2logL, logL,
			variances,
			propNugget=propNugget,			
			beta=betaHat, 
			sebeta,
			logDetVar=as.numeric(logDetVar),
			rpr=rpr)
	
	attributes(result)$varbetahat = varbetahat
	
	result

}

loglikGmrfOneRange = function(
		oneminusar=NULL, rangeInCells=NULL,
		Yvec, Xmat, NN, propNugget=0,
		shape=1,
		adjustEdges=FALSE,adjustParam=FALSE,
		adjustShape=FALSE,adjustMarginalVariance=FALSE) {
	
	Nx=attributes(NN)$Nx;Ny=attributes(NN)$Ny
	
	
	if(is.null(oneminusar)){

		NN = maternGmrfPrec(NN,
				param=c(shape=shape,
						range=as.vector(rangeInCells[1])),
				adjustEdges=adjustEdges,adjustParam=adjustParam,
				adjustShape=adjustShape,
				adjustMarginalVariance=adjustMarginalVariance)
		oneminusar=NA
	} else {
 		NN =  maternGmrfPrec(NN,
				param=c(shape=as.vector(shape),
						oneminusar=as.vector(oneminusar[1])),
				adjustEdges=adjustEdges,adjustParam=adjustParam,
				adjustShape=adjustShape,
				adjustMarginalVariance=adjustMarginalVariance)
		if(adjustParam) {
				rangeInCells = 	
					attributes(NN)$param$optimal['rangeInCells']
				shape=attributes(NN)$param$optimal['shape']
		} else {
			rangeInCells = 
					attributes(NN)$param$theo['rangeInCells']
		}
	}		

	cholPrec = Cholesky(NN,LDL=FALSE)
	theDet = 2*as.numeric(determinant(cholPrec,
					logarithm=TRUE)$modulus)
	
	Xp = solve(cholPrec,Xmat,system="P")
	Yp = solve(cholPrec,Yvec,system="P")

	
	LofQ = expand(cholPrec)$L
	XL = Matrix::crossprod(LofQ,Xp)
	YL = as.vector(
			Matrix::crossprod(LofQ,Yp)
	)

	argList = list(Yp=Yp,YL=YL,
			Xp=Xp,XL=XL,
			Qchol=cholPrec,
			detQ=theDet)


	if(!all(propNugget==0)) {
		argList$QLL = forceSymmetric(crossprod(LofQ,LofQ))
		argList$cholQLL = Cholesky(argList$QLL,LDL=FALSE)
		argList$YL = solve(argList$cholQLL, argList$YL,
				system='P')	
		argList$XL = solve(argList$cholQLL, argList$XL,
				system='P')	
		
	}
	
	#QLL=argList$QLL;cholQLL=argList$cholQLL;XL=argList$XL;YL=argList$YL
	
	# Q=NN;Qchol=cholPrec;detQ=theDet


# variance of Q inverse
midcellCoord = c(round(Nx*.5),round(Ny*0.5)) # the middle cell
midcell = c(Nx*(Ny-midcellCoord[2]) + midcellCoord[1])
midVec = sparseMatrix(midcell,1,x=1,
		dims=c(ncol(NN),1))

Qcentre = solve(cholPrec,
		solve(cholPrec, 
				solve(cholPrec,midVec,system="P"),
				system="A")	,
		system='Pt')[midcell]


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
	res$profL$oneminusar = cbind(
			oneminusar = x['oneminusar',1,],
			logL.ml=apply(x['logL.ml',,],2,max),
			logL.reml=apply(x['logL.reml',,],2,max)
	)
	res$profL$rangeInCells = list(
			ml = cbind(
				rangeInCells=	x['rangeInCells',1, ],
				logL=	apply(x['logL.ml',,],2,max)
			),
			reml = cbind(
					rangeInCells=x['rangeInCells',, ],
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
	
xsub = x[,grep("^se\\.",colnames(x),invert=TRUE)]
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
 
