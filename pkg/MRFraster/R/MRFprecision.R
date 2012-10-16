MRFmatrix = function(Nrow, Ncol, torus=FALSE){


	Ncells = Nrow * Ncol
	
	if(F) {
		mymat=matrix(1:Ncells, Nrow, Ncol,byrow=T)
		
	}
	
	
	
	library(spam)
	
	NNmatrix = spam(0, nrow=Ncells, ncol=Ncells)
# diagonals have 1
	diag(NNmatrix)=1


	# square left
	for(Drow in 1:Nrow) {
		for(Dcol in 2:Ncol){			
			thisCell = Dcol + Ncol*(Drow-1)
			NNmatrix[thisCell,thisCell-1] =
			NNmatrix[thisCell-1,thisCell] = 2
		}
	}

	# square up
	for(Drow in 2:Nrow) {
		for(Dcol in 1:Ncol){			
			thisCell = Dcol + Ncol*(Drow-1)
			NNmatrix[thisCell,thisCell-Ncol] =
					NNmatrix[thisCell-Ncol,thisCell] = 2
		}
	}

	# two up
	for(Drow in 3:Nrow) {
		for(Dcol in 1:Ncol){			
			thisCell = Dcol + Ncol*(Drow-1)
			NNmatrix[thisCell,thisCell-2*Ncol] =
					NNmatrix[thisCell-2*Ncol,thisCell] = 3
		}
	}
	
	
	# two left
	for(Drow in 1:Nrow) {
		for(Dcol in 3:Ncol){			
			thisCell = Dcol + Ncol*(Drow-1)
			NNmatrix[thisCell,thisCell-2] =
					NNmatrix[thisCell-2,thisCell] = 3
		}
	}

	# diag up left
	for(Drow in 2:Nrow) {
		for(Dcol in 2:Ncol){			
			thisCell = Dcol + Ncol*(Drow-1)
			newCell = thisCell - Ncol - 1
			NNmatrix[thisCell,newCell] =
					NNmatrix[newCell,thisCell] = 4
		}
	}

	# diag down left
	for(Drow in 2:Nrow) {
		for(Dcol in 1:(Ncol-1)){			
			thisCell = Dcol + Ncol*(Drow-1)
			newCell = thisCell - Ncol + 1
			NNmatrix[thisCell,newCell] =
					NNmatrix[newCell,thisCell] = 4
		}
	}

	# three to left
	for(Drow in 1:Nrow){
		for(Dcol in 4:Ncol){
			thisCell = Dcol + Ncol*(Drow-1)
			newCell = thisCell - 3  
			NNmatrix[thisCell, newCell]=
			NNmatrix[ newCell,thisCell]=5
			
		}
	}
	
	# three up
	for(Drow in 4:Nrow){
		for(Dcol in 1:Ncol){
			thisCell = Dcol + Ncol*(Drow-1)
			newCell = thisCell - 3*Ncol  
			NNmatrix[thisCell, newCell]=5
			NNmatrix[newCell,thisCell]=5
			
		}
	}
	
	# two down one left
	Cseq = 2:(Ncol)
	for(Drow in 1:(Nrow-2) ) {
		for(Dcol in Cseq){			
			thisCell = Dcol + Ncol*(Drow-1)
			newCell = thisCell + 2*Ncol - 1
			NNmatrix[thisCell,newCell] =
					NNmatrix[newCell,thisCell] = 6
		}
	}
	
	
	# two up one left
	Cseq = 2:(Ncol)
	for(Drow in 3:(Nrow) ) {
		for(Dcol in Cseq){			
			thisCell = Dcol + Ncol*(Drow-1)
			newCell = thisCell - 2*Ncol - 1
			NNmatrix[thisCell,newCell] =
					NNmatrix[newCell,thisCell] = 6
		}
	}
	
	
	# two left one up
	Cseq = 3:(Ncol)
	for(Drow in 2:(Nrow) ) {
		for(Dcol in Cseq){			
			thisCell = Dcol + Ncol*(Drow-1)
			newCell = thisCell - Ncol - 2
			NNmatrix[thisCell,newCell] =
					NNmatrix[newCell,thisCell] = 6
		}
	}
	
	
	# two left one down
Cseq = 3:(Ncol)
for(Drow in 1:(Nrow-1) ) {
	for(Dcol in Cseq){			
		thisCell = Dcol + Ncol*(Drow-1)
		newCell = thisCell + Ncol - 2
		NNmatrix[thisCell,newCell] =
				NNmatrix[newCell,thisCell] = 6
	}
}


if(F) {
	matrix(NNmatrix[,1], byrow=T, ncol=Ncol)
	matrix(NNmatrix[,12], byrow=T, ncol=Ncol)
	matrix(NNmatrix[,24], byrow=T, ncol=Ncol)
}

if(torus){
	
	# square left
	Dcol = 1
	for(Drow in 1:Nrow) {
			thisCell = Dcol + Ncol*(Drow-1)
			newCell = Ncells - Drow + 1
			NNmatrix[thisCell,thisCell-1] =
					NNmatrix[thisCell-1,thisCell] = 2
		
	}
	
	# square up
	for(Drow in 2:Nrow) {
		for(Dcol in 1:Ncol){			
			thisCell = Dcol + Ncol*(Drow-1)
			NNmatrix[thisCell,thisCell-Ncol] =
					NNmatrix[thisCell-Ncol,thisCell] = 2
		}
	}
	
	# two up
	for(Drow in 3:Nrow) {
		for(Dcol in 1:Ncol){			
			thisCell = Dcol + Ncol*(Drow-1)
			NNmatrix[thisCell,thisCell-2*Ncol] =
					NNmatrix[thisCell-2*Ncol,thisCell] = 3
		}
	}
	
	
	# two left
	for(Drow in 1:Nrow) {
		for(Dcol in 3:Ncol){			
			thisCell = Dcol + Ncol*(Drow-1)
			NNmatrix[thisCell,thisCell-2] =
					NNmatrix[thisCell-2,thisCell] = 3
		}
	}
	
	# diag up left
	for(Drow in 2:Nrow) {
		for(Dcol in 2:Ncol){			
			thisCell = Dcol + Ncol*(Drow-1)
			newCell = thisCell - Ncol - 1
			NNmatrix[thisCell,newCell] =
					NNmatrix[newCell,thisCell] = 4
		}
	}
	
	# diag down left
	for(Drow in 2:Nrow) {
		for(Dcol in 1:(Ncol-1)){			
			thisCell = Dcol + Ncol*(Drow-1)
			newCell = thisCell - Ncol + 1
			NNmatrix[thisCell,newCell] =
					NNmatrix[newCell,thisCell] = 4
		}
	}
	
	# three to left
	for(Drow in 1:Nrow){
		for(Dcol in 4:Ncol){
			thisCell = Dcol + Ncol*(Drow-1)
			newCell = thisCell - 3  
			NNmatrix[thisCell, newCell]=
					NNmatrix[ newCell,thisCell]=5
			
		}
	}
	
	# three up
	for(Drow in 4:Nrow){
		for(Dcol in 1:Ncol){
			thisCell = Dcol + Ncol*(Drow-1)
			newCell = thisCell - 3*Ncol  
			NNmatrix[thisCell, newCell]=5
			NNmatrix[newCell,thisCell]=5
			
		}
	}
	
	# two down one left
	Cseq = 2:(Ncol)
	for(Drow in 1:(Nrow-2) ) {
		for(Dcol in Cseq){			
			thisCell = Dcol + Ncol*(Drow-1)
			newCell = thisCell + 2*Ncol - 1
			NNmatrix[thisCell,newCell] =
					NNmatrix[newCell,thisCell] = 6
		}
	}
	
	
	# two up one left
	Cseq = 2:(Ncol)
	for(Drow in 3:(Nrow) ) {
		for(Dcol in Cseq){			
			thisCell = Dcol + Ncol*(Drow-1)
			newCell = thisCell - 2*Ncol - 1
			NNmatrix[thisCell,newCell] =
					NNmatrix[newCell,thisCell] = 6
		}
	}
	
	
	# two left one up
	Cseq = 3:(Ncol)
	for(Drow in 2:(Nrow) ) {
		for(Dcol in Cseq){			
			thisCell = Dcol + Ncol*(Drow-1)
			newCell = thisCell - Ncol - 2
			NNmatrix[thisCell,newCell] =
					NNmatrix[newCell,thisCell] = 6
		}
	}
	
	
	# two left one down
	Cseq = 3:(Ncol)
	for(Drow in 1:(Nrow-1) ) {
		for(Dcol in Cseq){			
			thisCell = Dcol + Ncol*(Drow-1)
			newCell = thisCell + Ncol - 2
			NNmatrix[thisCell,newCell] =
					NNmatrix[newCell,thisCell] = 6
		}
	}
	
	
}
	
	return(NNmatrix)
	
}
	
MRFprecision = function(sigmasq = 3,
		rho = 2,
		nu=1,
		gridSize = 1,
		Nrow=30, Ncol=30, 
		returnCov = F) {

	
	NNmatrix = MRFmatrix(Nrow, Ncol)
	
	kappa = sqrt(8*nu)/rho
	a =  (gridSize*kappa)^2 + 4  
	asq = a*a
	
	if(nu==1){
	mrfCoef = c(diag=4 + asq, 
		firstN=-2*a, 
		twoNcorner=2, twoNside=1)
	 mrfCoef = mrfCoef /  (sigmasq*4*pi*nu*(a-4)^nu)
	
	mrfCoefSub = c(
			mrfCoef[c("diag","firstN","twoNcorner","twoNside")], 
			rep(0, 2))
} else if(nu==2){
	
	mrfCoef = c(
			diag=a*(asq+12),
			first=-3*(asq+3),
			twoCorner=6*a,
			twoSide=3*a,
			threeSide=-1,
			threeCorner=-3
			)
	#		mrfCoef = mrfCoef /  (sigmasq*4*pi*nu*(a-4)^nu)
	mrfCoefSub = mrfCoef[c(
					"diag","first","twoCorner","twoSide",
					"threeSide","threeCorner")]		
	
} else {
warning("MRF precision matrix only defined for nu = 1 or 2")	
}
	
	precMat = NNmatrix
	
	precMat@entries = mrfCoefSub[NNmatrix@entries]
	
 
	
if(returnCov) {
	result = list(prec=precMat)
	
	xseq= seq(0,by=gridSize, length=Ncol)
	yseq = seq(to=0,by=-gridSize, length=Nrow)
	
	gridMat = outer(yseq*1i, xseq, FUN="+")
	
	
	gridC = as.vector(t(gridMat))
	
	distmat = Mod(outer(gridC, gridC, FUN="-"))*kappa
	
	covMat = distmat^nu *besselK(distmat, nu)
	diag(covMat)=1
	covMat = covMat * (sigmasq/(2^(nu-1)*gamma(nu)) )
	
	result$cov = covMat
	
} else {
result=precMat
}


if(F) {
	
	temp= precMat  %*% result$cov
	
	D=15
	matrix(c(precMat[D,]), byrow=T, ncol=Ncol)
	
	matrix(c(precMat[1,]), byrow=T, ncol=Ncol)
	matrix(c(NNmatrix[1,]), byrow=T, ncol=Ncol)
	
	matrix(c(precMat[6,]), byrow=T, ncol=Ncol)
	
	matrix(c(precMat[30,]), byrow=T, ncol=Ncol)
	
	matrix(c(precMat[25,]), byrow=T, ncol=Ncol)
	
	matrix(solve(covMat)[25,], byrow=T, ncol=Ncol)
	
	
	for(D in 1:10) {
		cat(D)
		print(matrix(c(NNmatrix[D,]), byrow=T, ncol=Ncol))
		scan()
	}
	
}

result
	
}