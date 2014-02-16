conditionalGmrf = function(param, MoreArgs=NULL,
		Yvec=MoreArgs$Yvec,Xmat=MoreArgs$Xmat,
		NN=MoreArgs$NN,mc.cores=1) {
	
	ar = param['ar']
	maternShape=param['maternShape']
	tausq = param['tausq']
	xisq = param['xisq']
	
	
	
	if(maternShape==0) {
		paramXX = c(1,-ar/4, rep(0,4))
		NN@x = paramXX[NN@x]
	} else if(maternShape==1){
		paramXX = c(4 + (4/ar)^2, -8/ar, 2, 1, rep(0,2))
		NN@x = paramXX[NN@x]
	} else if(maternShape==2){
		paramXX = c((4/ar)*(12+16/ar^2),-3*(3 +  16*ar^2), 24/ar,12/ar, -3, -1)
		NN@x = paramXX[NN@x]
	}else {
		warning("shape should be 0, 1, or 2")
	}
	
	IcQ=NN
	IcQ@x = IcQ@x * (tausq/xisq)
	diag(IcQ) = diag(IcQ) +1 
	
	fixed=as.vector(Xmat %*% 
			param[paste("beta",colnames(Xmat), sep=".")]
)
	resids = Yvec -	fixed

	
	if(mc.cores==1) {
		cIcQ = Cholesky(IcQ)
		cQ = Cholesky(NN)
	} else {
		cIcQ = mcparallel(Cholesky(IcQ),name='cIcQ')
		cQ = mcparallel(Cholesky(NN),name='cQ')
		thec = mccollect(list(cIcQ, cQ))
		cIcQ = thec[['cIcQ']]
		cQ = thec[['cQ']]
	}

	EUY = as.vector(solve(cIcQ, resids))
	
	varOneCell = function(D) {
		thisD = sparseMatrix(D,1,x=1,dims=c(Ny,1))
		solveQ = solve(cQ, thisD)
		as.vector(
				solveQ[D] - sum(solve(cIcQ, thisD) * solveQ)
		)
	}
	newenv =new.env()
	junk=eval(library("Matrix"),envir=newenv)
	assign("Ny",length(resids),envir=newenv)
	assign("cIcQ",cIcQ,envir=newenv)
	assign("cQ",cQ,envir=newenv)
	environment(varOneCell) = newenv
	

	if(mc.cores==1) {
		thediag = mapply(varOneCell, 1:length(resids))
	} else {
		thediag = mcmapply(varOneCell, 1:length(resids),
				mc.cores=mc.cores,SIMPLIFY=TRUE)
	}

	VUY = xisq * thediag

	result=cbind(random=EUY,krigeSd=sqrt(VUY),
			fixed,fitted=fixed+EUY,resids)
	result
}

 