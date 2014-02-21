conditionalGmrf = function(param,
		Yvec,Xmat,
		predRaster, NN=NNmat(predRaster),
		mc.cores=1, ...) {
	
	ar = 1-param['oneminusar']
	maternShape=param['shape']
	tausq = param['tausq.ml']
	sigsq = param['sigmasq.ml']
	
	
	fixed=as.vector(Xmat %*% 
					param[paste("beta",colnames(Xmat), sep=".")]
	)
	
	Q = maternGmrfPrec(NN,param,...)
	Qchol = Cholesky(Q, LDL=FALSE)

	residsOrig = Yvec -	fixed
	resids = solve(Qchol,
			residsOrig,
			system='P')
	
	nuggetSigsq = tausq/sigsq
	# nuggetSigsq = tausq/sigsq

	LofQ = expand(Qchol)$L
	QLL = forceSymmetric(crossprod(LofQ,LofQ))
	diag(QLL) = diag(QLL) + 1/nuggetSigsq
	cholIcQ = Cholesky(QLL, LDL=FALSE)

	residsP = solve(cIcQ, resids,system='P')
	
	
	EUY = solve(Qchol,
			solve(cIcQ,
			as.vector(solve(cIcQ, residsP,system='A')),
			system="Pt"),
	system='Pt')
	
	varOneCell = function(D) {
		thisD = sparseMatrix(D,1,x=1,dims=c(Ny,1))
		solveQ = solve(Qchol, thisD,system='A')
		solveQp = solve(cIcQ, solveQ,system='P')
		thisDp = solve(cIcQ, thisDp ,system='P')
		as.vector(
				solveQ[D] - sum(solve(cIcQ, thisDp,
								system='A') * solveQ)
		)
	}

	
	if(mc.cores==1) {
		thediag = mapply(varOneCell, 1:length(resids))
	} else {
		thediag = mcmapply(varOneCell, 1:length(resids),
				mc.cores=mc.cores,SIMPLIFY=TRUE)
	}

	VUY = sigsq * solve(
			cIcQ,
			thediag,
			system='Pt'
	)

	result=cbind(random=EUY,krigeSd=sqrt(VUY),
			fixed,fitted=fixed+EUY,residsOrig)
	result
}

 