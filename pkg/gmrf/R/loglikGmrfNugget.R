logLikGmrfNugget = function(ar, propNugget=5, Yvec, Xmat, 
		NN, maternShape=1) {
	

	if(length(propNugget)==1) {
	if(round(propNugget)==propNugget & propNugget>1 & propNugget>0){
		# create a sequence
		propNugget = seq(0,1,len=propNugget+1)[
				-propNugget-1
		]
	}
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
	
	if(length(nuggetDivXisq)==1) {
		argList$nuggetDivXisq=nuggetDivXisq
		res=do.call(
				logLikGmrfNuggetSingle,
				argList)
	} else {
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
	}
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

XprecXinv = solve(Matrix::crossprod(Xbreve,Xbreve)) 

betaHat = as.vector(
		XprecXinv %*% 
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


result =c(m2logL=m2logL, logL=NA,beta=betaHat, 
		sigmasq=NA, rangeInCells=NA,
		tausq=nuggetVarHat,ar=NA,
		propNugget=NA, maternShape=NA,
		xisq=xisqHat)

varbetahat = as.matrix(XprecXinv*xisqHat)
dimnames(varbetahat) = list(names(betaHat),names(betaHat))
attributes(result)$varbetahat = varbetahat


result

}





summaryGmrfFitNugget = function(applyResult,MoreArgs,
		fun=logLikGmrfNugget) {
	
	if(is.list(applyResult))
		applyResult = simplify2array(applyResult)
	
	x = aperm(applyResult,3:1)
	x = matrix(x, ncol=dim(x)[3])
	colnames(x) = dimnames(applyResult)[[1]]
	
	mleRow = which.max(x[,'logL'])
	mle = x[mleRow,]
	

	MoreArgs$ar = mle['ar']
	MoreArgs$propNugget = mle['propNugget']
	
	varBetaHat = attributes(
			do.call(fun, MoreArgs) 
		)$varbetahat

		
	thebetas = grep("^beta",names(mle), value=TRUE)
		
	betaMLE = mle[thebetas]

	
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
		withinCI = which(x[,'m2logL'] - 
						mle['m2logL'] < qchisq(0.95, 1))
		
		parCI= t(apply(
						x[withinCI, ,drop=FALSE],
						2,range))
		colnames(parCI) = c('q0.025', 'q0.975')

		
		
		parMat = cbind(mle=mle,parCI,se=NA,pval=NA)
		parMat[rownames(betamat),colnames(betamat)] = betamat
		
		
		
		result = list(
						summary=parMat,
						profL = list(
								)
				)

		for(D in c('rangeInCells','ar','propNugget')) {
			theunique = sort(unique(x[,D]))
			thefac = list(factor(x[,D],levels=theunique))
			
			result$profL[[D]] = 
					aggregate(x[,'logL'], thefac, min)
			colnames(result$profL[[D]]) = c('x','y')
		}
				
		result
}