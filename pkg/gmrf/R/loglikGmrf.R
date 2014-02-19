loglikGmrfSingle = function(
		propNugget, Yvec, YL, 
		Xmat, XL,Q, Qchol, detQ) {
	
	if(propNugget>0){
		
		# propNugget = tausq / (tausq + sigsq)
		nuggetSigsq = 1/(1/propNugget -1)
		
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
		
		resids =  Ybreve -  Xbreve %*% betaHat 
		
		sigsq = as.numeric(
				Matrix::crossprod(resids)/length(Ytilde)
		)
		
		nugget = sigsq * nuggetSigsq
		
		m2logL = as.numeric(
				-2*determinant(cholIcQ,logarithm=TRUE)$modulus -
						detQ + length(Yvec)*log(sigsq) 
		)
		varbetahat=NULL
	} else { # no nugget 
		
		XprecX = Matrix::tcrossprod(XL,XL)
		XprecXinv = solve(XprecX)
		
		betaHat = as.numeric(XprecXinv %*% 
						Matrix::tcrossprod(XL , YL))
		names(betaHat) = colnames(Xmat)
		
		# resids ~ N(0,sigsq*prec^(-1))
		resids = as.vector(Yvec - Xmat %*% betaHat)
		residsL = Matrix::crossprod(
				solve(Qchol,resids,system="P"), 
				expand(Qchol)$L)
		rpr = as.numeric(Matrix::tcrossprod(residsL,residsL))
		sigsq = rpr/ length(Yvec)
		
		varbetahat = as.matrix(XprecXinv*sigsq)
		dimnames(varbetahat) = list(names(betaHat),names(betaHat))
		
		m2logL = length(Yvec)*log(sigsq) - detQ
		nugget=0
	}
	
	
	
	
	result =c(m2logL=m2logL, logL=-m2logL/2,
			beta=betaHat, 
			sigmasq=sigsq, 
			tausq=nugget, 
			propNugget=propNugget)
	
	attributes(result)$varbetahat = varbetahat
	
	
	result
	
}

loglikGmrf = function(oneminusar=NULL, rangeInCells=NULL,
		Yvec, Xmat, NN, propNugget=0,
		shape=1,
		adjustEdges=FALSE,adjustParam=FALSE) {
	
	if(is.null(oneminusar)){

		NN = geostatsp::maternGmrfPrec(NN,
				param=c(shape=shape,
						range=oneminusar),
				adjustEdges=adjustEdges,adjustParam=adjustParam)
	} else {
		aroveronemar = (1-oneminusar)/oneminusar
		NN = geostatsp::maternGmrfPrec(NN,
				param=c(variance=1,shape=maternShape,cellSize=1,
						onemar=as.vector(1-ar)),
				adjustEdges=adjustEdges,adjustParam=adjustParam)
		
		if(adjustParam) {
			rangeInCells = 	attributes(NN)$param$sameShape['rangeInCells']
		} else {
			rangeInCells = attributes(NN)$param$theo['rangeInCells']
		}
	}

	if(!any(maternShape==c(0,1,2))) {
		warning("shape should be 0, 1, or 2")
	}
	
	
	
	cholPrec = Cholesky(NN,LDL=FALSE)
	theDet = 2*as.numeric(determinant(cholPrec,
					logarithm=TRUE)$modulus)
	
	# beta hat = (X'prec/C X)^{-1} X'prec/C y
	LofQ = expand(cholPrec)$L
	XL = Matrix::crossprod(solve(cholPrec,Xmat,system="P"),
			LofQ)
	YL = Matrix::crossprod(
			solve(cholPrec,Yvec,system="P"),LofQ)
	
	if(length(propNugget)==1) {
		res = loglikGmrfSingle(
				propNugget, Yvec, YL, 
				Xmat, XL,Q, cholPrec, theDet) 	
		res = c(res,
				rangeInCells=
						as.numeric(rangeInCells),
				ar=as.numeric(ar),
				maternShape=as.numeric(maternShape)
			)
		
	} else {
		

	argList = list(Yvec=Yvec,YL=YL,
			Xmat=Xmat,XL=XL,
			Q=NN,Qchol=cholPrec,
			detQ=theDet)

	res = mapply(loglikGmrfSingle,
			propNugget=propNugget,			
			MoreArgs=argList,
			SIMPLIFY=TRUE
	)
	res = rbind(res,
		rangeInCells=as.numeric(rangeInCells),
		ar=as.numeric(ar),
		maternShape=as.numeric(maternShape)
	)
	
}

res 

}

bob=function(){
	
	XprecX = Matrix::tcrossprod(XL,XL)
	XprecXinv = solve(XprecX)

	betahat = as.numeric(XprecXinv %*% 
					Matrix::tcrossprod(XL , YL))
	names(betahat) = colnames(Xmat)
	
	
	# resids ~ N(0,sigsq*prec^(-1))
	resids = as.vector(Yvec - Xmat %*% betahat)
	residsL = Matrix::crossprod(
			solve(cholPrec,resids,system="P"), 
			cholPrec)
	rpr = as.numeric(Matrix::tcrossprod(residsL,residsL))
	sigsq = rpr/ length(Yvec)

	varbetahat = as.matrix(XprecXinv*sigsq)
	dimnames(varbetahat) = list(names(betahat),names(betahat))
	
	m2logL = length(Yvec)*log(sigsq) - theDet
	
	result = c(
			m2logL = m2logL,
			logL = -m2logL/2,
			rangeInCells=rangeInCells,
			sigmasq= sigsq,
			beta=betahat,
			maternShape=maternShape,
			ar=ar,
			logdet=theDet
	)
	
	attributes(result)$varbetahat = varbetahat

	result
}
 
summaryGmrfFit = function(applyResult,MoreArgs) {
	
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
 
