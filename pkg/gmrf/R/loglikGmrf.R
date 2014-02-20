loglikGmrfGivenQ = function(
		propNugget, Yvec, YL, 
		Xmat, XL,Q, Qchol, detQ,
		Qcentre) {
	
	if(propNugget>0){
		
		# propNugget = tausq / (tausq + sigsq)
		nuggetSigsq = 1/( (1/propNugget) -1)
		
		IcQ = Q
		IcQ@x = IcQ@x * nuggetSigsq 
		diag(IcQ) = diag(IcQ) +1 
		cholIcQ = Cholesky(IcQ,LDL=FALSE)
		
		Ybreve = solve(cholIcQ, 
				solve(cholIcQ,YL,system="P"),
				system="L")
		Xbreve = solve(cholIcQ, 
				solve(cholIcQ,XL,system="P"),
				system="L")	
		
		XprecXinv = solve(Matrix::crossprod(Xbreve,Xbreve)) 
		
		betaHat = as.vector(
				XprecXinv %*% 
						Matrix::crossprod(Xbreve,Ybreve) 
		)
		names(betaHat) = colnames(Xmat)
		
		residsL =  Ybreve -  Xbreve %*% betaHat 
		
		
		thedet = as.numeric(
				-2*determinant(cholIcQ,logarithm=TRUE)$modulus -
						detQ 
		)
		
		varbetahat=NULL
		
	
	} else { # no nugget 
		nuggetSigsq = 0
		
		XprecX = Matrix::crossprod(XL,XL)
		XprecXinv = solve(XprecX)
		
		betaHat = as.numeric(XprecXinv %*% 
						Matrix::crossprod(XL , YL))
		names(betaHat) = colnames(Xmat)
		
		# resids ~ N(0,sigsq*prec^(-1))
		resids = as.vector(Yvec - as.vector(Xmat %*% betaHat))
		residsL = Matrix::crossprod(
				solve(Qchol,resids,system="P"), 
				expand(Qchol)$L)
		
		logDetVar=-detQ
		# matrix operations
		# Qchol, cholCovMat = Matrix::chol(covMat)
		# XL, cholCovInvX = Matrix::solve(cholCovMat, covariates)
		# XprecX, cholCovInvXcross = Matrix::crossprod(cholCovInvX)
		# XprecXinv, cholCovInvXcrossInv = Matrix::solve(cholCovInvXcross)
	}

	rpr = as.numeric(Matrix::tcrossprod(residsL,residsL))

	N=c(ml=0,reml=-ncol(Xmat))+length(Yvec)
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
		rep(Qcentre*rpr/N,c(2,2)) 
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
			along=3
			)
	dimnames(varbetahat) = 
					list(names(betaHat),names(betaHat),
							names(constHat))

			
	result =c(m2logL, logL,
			variances,
			propNugget=propNugget,			
			beta=betaHat, 
			sebeta)
	
	attributes(result)$varbetahat = varbetahat
	
	result

}

loglikGmrfOneRange = function(
		oneminusar=NULL, rangeInCells=NULL,
		Yvec, Xmat, NN, propNugget=0,
		shape=1,
		adjustEdges=FALSE,adjustParam=FALSE,
		adjustShape=FALSE) {
	
	Nx=attributes(NN)$Nx;Ny=attributes(NN)$Ny
	
	
	if(is.null(oneminusar)){

		NN = maternGmrfPrec(NN,
				param=c(shape=shape,
						range=as.vector(rangeInCells[1])),
				adjustEdges=adjustEdges,adjustParam=adjustParam,
				adjustShape=adjustShape)
		oneminusar=NA
	} else {
 		NN =  maternGmrfPrec(NN,
				param=c(shape=shape,
						oneminusar=oneminusar[1]),
				adjustEdges=adjustEdges,adjustParam=adjustParam,
				adjustShape=adjustShape)
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
	
	# beta hat = (X'prec/C X)^{-1} X'prec/C y
	LofQ = expand(cholPrec)$L
	XL = t(Matrix::crossprod(solve(cholPrec,Xmat,system="P"),
			LofQ))
	YL = as.vector(Matrix::crossprod(
			solve(cholPrec,Yvec,system="P"),LofQ)
	)	

	# variance of Q inverse
	midcellCoord = c(round(Nx*.5),round(Ny*0.5)) # the middle cell
	midcell = c(Nx*(Ny-midcellCoord[2]) + midcellCoord[1])
	midVec = sparseMatrix(midcell,1,x=1,
			dims=c(ncol(NN),1))
	
	Qcentre = max(
			solve(cholPrec, 
			solve(cholPrec,midVec,system="P"),
			system="A")	
	)
	
	argList = list(Yvec=Yvec,YL=YL,
			Xmat=Xmat,XL=XL,
			Q=NN,Qchol=cholPrec,
			detQ=theDet,Qcentre=Qcentre)
	
	if(length(propNugget)==1) {
		
		argList$propNugget=propNugget

		res = do.call(loglikGmrfGivenQ,argList)

		res = c(res,
				rangeInCells=
						as.numeric(rangeInCells),
				oneminusar=as.numeric(oneminusar),
				shape=as.numeric(shape)
			)
		
	} else {
		

		argList = list(Yvec=Yvec,YL=YL,
			Xmat=Xmat,XL=XL,
			Q=NN,Qchol=cholPrec,
			detQ=theDet,Qcentre=Qcentre)

		res = mapply(loglikGmrfGivenQ,
			propNugget=propNugget,			
			MoreArgs=argList,SIMPLIFY=TRUE
		)
		res = rbind(res,
			rangeInCells=as.numeric(rangeInCells),
			oneminusar=as.numeric(oneminusar),
			shape=as.numeric(shape)
		)
	
}

res

}



loglikGmrf = function(
		oneminusar=NULL, rangeInCells=NULL,
		Yvec, Xmat, NN, propNugget=0,
		shape=1,
		adjustEdges=FALSE,adjustParam=FALSE,
		adjustShape=FALSE,
		mc.cores=1) {

	
		argList = list(Yvec=Yvec, 
			Xmat=Xmat, NN=NN,
			propNugget=propNugget,
			shape=shape,
			adjustEdges=adjustEdges,
			adjustParam=adjustParam,
			adjustShape=adjustShape
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
			res= simplify2array(res)
		}
		res
}

 


summaryGmrfFit= function(applyResult,MoreArgs) {
	UseMethod("summaryGmrfFit")	
}
 
summaryGmrfFit.array = function() {
	
}

summaryGmrfFit.matrix = function(applyResult,MoreArgs) {
	
	fun=loglikGmrf
	
	applyResult = t(applyResult)

	MLErow = applyResult[which.max(applyResult[ ,'logL']),]	
	
	varBetaHat = attributes(
			do.call(fun,c(list(ar=MLErow["ar"]),MoreArgs))
	)$varbetahat
	
	thebetas = grep("^beta",names(MLErow), value=TRUE)
		
	betaMLE = MLErow[thebetas]	
	thebetas2 = gsub("^beta\\.","",thebetas)
	
	betase = sqrt(diag(varBetaHat[thebetas2,thebetas2]))
	
	betamat = cbind(mle=betaMLE,
			se=betase,
			q0.025=betaMLE - 2*betase,
			q0.975=betaMLE + 2*betase,
			pval = 2*pnorm(
					abs(betaMLE)/
							betase,
					lower=FALSE)	
	)
	withinCI = which(applyResult[,'m2logL'] - 
					MLErow['m2logL'] < qchisq(0.95, 1))
	
	parCI= t(apply(
					applyResult[withinCI, ,drop=FALSE],
					2,range))
	colnames(parCI) = c('q0.025', 'q0.975')
	
	parMat = cbind(mle=MLErow,parCI,se=NA,pval=NA)
	parMat[rownames(betamat),colnames(betamat)] = betamat

	return(list(
					summary=parMat,
					profL = cbind(
							applyResult, 
							logL = -applyResult[,'m2logL']/2)
			)
	)		
}
 
