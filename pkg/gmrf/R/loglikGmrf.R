logLikGmrf = function(ar, Yvec, Xmat, NN, maternShape=0) {
	

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
	# beta hat = (X'prec/C X)^{-1} X'prec/C y
	Xprec = Matrix::crossprod(Xmat,NN)
	XprecX = Xprec %*% Xmat
	XprecXinv = solve(XprecX)
	betahat = as.numeric(XprecXinv %*% (Xprec %*% Yvec))
	resids = as.vector(Yvec - Xmat %*% betahat)
	
	# resids ~ N(0,C*prec^(-1))
	# pr(resids) = C^(-n/2) * sqrt(det(prec)) * exp(-resid prec resid / 2 C)
	# - 2 log pr(resids) = n log C - log(det(prec)) + resid prec resid /  C
	#  deriv = 0 when  n/ C = rpr C^(-2), Chat = rpr/n
	rpr = as.numeric(crossprod(resids,NN)%*%resids)
	Chat = rpr/ length(Yvec)
	
	#var betahat =  xpxinv x' prec  var(Y) prec x xpxinv'
	#var betahat =  xpxinv x'  prec x xpxinv' C
	#var betahat =  xpxinv' C
	theDet = as.numeric(determinant(NN,logarithm=TRUE)$modulus)
	Chat = as.numeric(Chat)
	
	names(betahat) = colnames(Xmat)
	
	result = c(
			m2logL = as.numeric(length(Yvec) * log(Chat) - theDet ),
			rangeInCells=sqrt(ar/(1-ar))*sqrt(2*maternShape),
			sigma=exp(
					log(Chat) - 
							(maternShape+1)*log(4)+
						maternShape*(log(ar)-log(1-ar)) - 
						log(pi)),
			beta=betahat,
			maternShape=maternShape,
			ar=ar,
			xi=Chat,
			logdet=theDet
	)
	
	varbetahat = as.matrix(XprecXinv*Chat)
	dimnames(varbetahat) = list(names(betahat),names(betahat))
	attributes(result)$varbetahat = varbetahat

	result
}
 
summaryGmrfFit = function(applyResult,MoreArgs,fun=logLikGmrf) {
	
	
	applyResult = t(applyResult)

	MLErow = applyResult[which.min(applyResult[ ,'m2logL']),]	
	
	varBetaHat = attributes(do.call(fun,c(list(ar=MLErow["ar"]),MoreArgs)))$varbetahat
	



	thebetas = grep("^beta",names(MLErow), value=TRUE)
		
	betaMLE = MLErow[thebetas]
	betase = sqrt(diag(varBetaHat[thebetas,thebetas]))
	
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
 
