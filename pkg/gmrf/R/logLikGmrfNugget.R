
#+ logLikGmrfNugget
logLikGmrfNugget = function(param, Yvec, Xmat, NN) {
	
	if(length(param)==1) {
		NN@x = c(1,-param)[NN@x]
	} else {
		NN@x = param[NN@x]
	}
	
	Q = NN
	cholQ = chol(Q)
	
	Ytilde = solve(cholQ, Yvec)
	Xtilde = solve(cholQ, Xmat)
	
	ItauQ = nugget * Q
	diag(ItauQ) = diag(ItauQ)+1
	cholItauQ = chol(ItauQ)
	
	Xstuff = solve(cholItauQ, Xtilde)
	Ystuff = solve(cholItauQ, Ytilde)
	
	XprecX = crossprod(Xstuff,Xstuff)
	betahat = solve(XprecX ) %*% crossprod(Xstuff,Ystuff)
	
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
	
}
#'