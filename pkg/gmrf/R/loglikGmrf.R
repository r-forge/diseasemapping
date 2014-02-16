loglikGmrf = function(ar=NULL, rangeInCells=NULL,
		Yvec, Xmat, NN, propNugget=0,
		maternShape=1,
		adjustEdges=FALSE) {
	
	if(is.null(ar)){
		aroveronemar = rangeInCells^2/(2*maternShape)
		ar=aroveronemar/(1+aroveronemar)
	} else {
		aroveronemar = ar/(1-ar)
		rangeInCells = sqrt(2*maternShape*aroveronemar)
	}

	if(!any(maternShape==c(0,1,2))) {
		warning("shape should be 0, 1, or 2")
	}
	
	
	NN = geostatsp::maternGmrfPrec(NN,
			param=c(variance=1,shape=maternShape,cellSize=1,
					range=rangeInCells),
			adjustEdges=adjustEdges)
	cholPrec = Cholesky(NN,LDL=FALSE)
	theDet = 2*as.numeric(determinant(cholPrec,
					logarithm=TRUE)$modulus)
	
	# beta hat = (X'prec/C X)^{-1} X'prec/C y
	XL = Matrix::crossprod(solve(cholPrec,Xmat,system="P"),
			cholPrec)
	YL = Matrix::crossprod(
			solve(cholPrec,Yvec,system="P"),cholPrec)
	
	
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
	
	varBetaHat = attributes(do.call(fun,c(list(ar=MLErow["ar"]),MoreArgs)))$varbetahat
	



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
 
