logLikGmrfNugget = function(ar, propNugget=5, Yvec, Xmat, 
		NN, maternShape=1) {
	
	if(round(propNugget)==propNugget & propNugget>1){
		propNugget = seq(0,1,len=propNugget+1)[
				-propNugget-1
		]
	}
	
	varSpatialDivXisq = exp(
			maternShape*(log(ar)-log(1-ar))-
						log(pi)-
						(maternShape+1)*log(4)
						)
	
	nuggetDivXisq = 1/( (1/propNugget)-1) * varSpatialDivXisq
	
	
	if(maternShape==0) {
		param = c(1,-ar/4, rep(0,4))
		NN@x = param[NN@x]
	} else if(maternShape==1){
		param = c(4 + (4/ar)^2, -8/ar, 2, 1, rep(0,2))
		NN@x = param[NN@x]
	} else if(maternShape==2){
		param = c((4/ar)*(12+16/ar^2),-3*(3 +  16*ar^2), 24/ar,12/ar, -3, -1)
		NN@x = param[NN@x]
	}else {
		warning("shape should be 0, 1, or 2")
	}

	cholQ =  Matrix::chol(NN)
	detCholQ = determinant(cholQ,logarithm=TRUE)$modulus
	
	Ytilde = Matrix::crossprod(cholQ, Yvec)
	Xtilde = Matrix::crossprod(cholQ, Xmat)

	argList = list(Ytilde=Ytilde,Xtilde=Xtilde,
			Q=NN,detCholQ = detCholQ)
	

	res = mapply(logLikGmrfNuggetSingle,
			nuggetDivXisq=nuggetDivXisq,			
			MoreArgs=argList,
			SIMPLIFY=TRUE
	)
	

	
	res["sigmasq",] = res["xisq",] * varSpatialDivXisq
	res["rangeInCells",]=sqrt(ar/(1-ar))*sqrt(2*maternShape)
	res["ar",]=ar
	res["propNugget",]=propNugget
	res["maternShape",]=maternShape
	res["logL",]=-res["m2logL",]/2
	
	res
	
}

logLikGmrfNuggetSingle = function(
		nuggetDivXisq, Ytilde, 
		Xtilde,Q, detCholQ) {
IcQ = nuggetDivXisq * Q
diag(IcQ) = diag(IcQ) +1 
cholIcQ = chol(IcQ)

Ybreve = solve(cholIcQ, Ytilde)	
Xbreve = solve(cholIcQ, Xtilde)	

betaHat = as.vector(
		solve(Matrix::crossprod(Xbreve,Xbreve)) %*% 
		Matrix::crossprod(Xbreve,Ybreve) 
)
names(betaHat) = colnames(Xtilde)

resids = Ybreve -  Xbreve %*% betaHat

xisqHat = as.numeric(
		Matrix::crossprod(resids)/length(Ytilde)
	)

m2logL = as.numeric(
		2*determinant(cholIcQ,logarithm=TRUE)$modulus -
		2* detCholQ + length(Ytilde)*log(xisqHat) 
	)

nuggetVarHat = xisqHat * nuggetDivXisq

c(m2logL=m2logL, logL=NA,beta=betaHat, 
		sigmasq=NA, rangeInCells=NA,
		tausq=nuggetVarHat,ar=NA,
		propNugget=NA, maternShape=NA,
		xisq=xisqHat)

}