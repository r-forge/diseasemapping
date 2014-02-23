conditionalGmrf = function(param,
		Yvec,Xmat, NN,
		template=NULL, 
		mc.cores=1, ...) {
	
	names(param) = gsub("sigmasq","variance",names(param))
	names(param) = gsub("tausq","nugget",names(param))
	

	fixed=as.vector(Xmat %*% 
					param[paste("beta",colnames(Xmat), sep=".")]
	)
	
	
	Q = maternGmrfPrec(NN,param[c('oneminusar','shape','variance')],...)
	Qchol = Cholesky(Q, LDL=FALSE)

	LofQ = expand(Qchol)$L
	pRev = as(ncol(LofQ):1, "pMatrix")		
	lQLL =  as( t(LofQ %*% pRev) %*% pRev,'dtCMatrix')
	QLL = tcrossprod(lQLL)
	diag(QLL) = diag(QLL) + param['nugget']
	QLL = forceSymmetric(QLL)
	
	cholIcQ = Cholesky(QLL,LDL=FALSE,perm=TRUE)
	
	# QLL = P' L  L' P 
	# QLLorig = Prev P' L   L' P Prev'
	
	ptwice =   as(expand(cholQLL)$P,'sparseMatrix') %*% 
			t(as(pRev,'sparseMatrix'))
	
	ptwice2 = as(ptwice, 'pMatrix')
	cholIcQ@perm = as.integer(ptwice2@perm-1)

	residsOrig = Yvec -	fixed
	
	varyinvresid = as.vector(solve(cholIcQ, residsOrig,system='A'))
	
	EUY = as.vector(solve(Qchol,
					varyinvresid))
	
theidm=c(length(EUY),1)
	varOneCell = function(D) {
		thisD = sparseMatrix(D,1,x=1,dims=theidm)
		solveQ = solve(Qchol, thisD)
#		solveQp = solve(cholIcQ, solveQ,system='P')
#		thisDp = solve(cholIcQ, thisD ,system='P')
		as.vector(
				solveQ[D] - sum(solve(cholIcQ, thisD,
								system='A') * solveQ)
		)
	}

	
	if(mc.cores==1) {
		thediag = mapply(varOneCell, 1:length(resids))
	} else {
		thediag = parallel::mcmapply(varOneCell, 1:length(EUY),
				mc.cores=mc.cores,SIMPLIFY=TRUE)
	}

	
	
	VUY = as.vector(
		param['variance'] * 
			thediag
	)

	
	result=cbind(random=EUY,krigeSd= sqrt(VUY),
			fixed=fixed,fitted=fixed+EUY,resids=residsOrig)
	if(!is.null(template)){
		resRast = raster::brick(raster(template), nl=ncol(result))
		names(resRast) = colnames(result)
		values(resRast) = as.vector(result)
		result = list(df=result,rast=resRast)		
		
	}
	result
}

 