library(geostatsp)
matrix(NNmat(7, 7)[,25], 7, 7)


params=c(range = 3,
		cellSize=0.5,
		shape=1,
		variance=50^2)

theNN = NNmat(40,20,nearest=params['shape']+1)
# precision matrix without adjusting for edge effects
precMat =maternGmrfPrec(theNN, param=params) 
# and with the adjustment
precMatCorr =maternGmrfPrec(theNN, param=params, 
		adjustEdges=TRUE) 

#
N=theNN;Ny=N;param=params

Nx = attributes(theNN)$Nx
Ny = attributes(theNN)$Ny
midcell = Nx*round(Ny/2) + round(Nx/2) # the middle cell
edgecell = Nx*5 + 5 # cell near corner

# show precision of cell 32,32 
precMid=matrix(precMat[,midcell], Ny, Nx, byrow=TRUE)
precMid[round(Ny/2)+seq(-3, 3), round(Nx/2)+seq(-3, +3)]

 
# variance matrices
Ncell = Nx*Ny
midVec = sparseMatrix(midcell,1,x=1,dims=c(Ncell,1))
edgeVec = sparseMatrix(edgecell,1,x=1,dims=c(Ncell,1))

therast = attributes(precMat)$raster
values(therast)=NA
varRast = brick(therast, therast,therast,therast)
names(varRast) = c('mid','edge','midCor','edgeCor')
values(varRast[['mid']]) = as.vector(Matrix::solve(precMat, midVec))
values(varRast[['edge']]) = as.vector(Matrix::solve(precMat, edgeVec))
values(varRast[['midCor']]) =  as.vector(Matrix::solve(precMatCorr, midVec))
values(varRast[['edgeCor']]) = as.vector(Matrix::solve(precMatCorr, edgeVec))

projection(varRast) = "+init:units=m"


if(!interactive()){
	pdf("maternGmrfPredRasters.pdf",height=8,width=5)
}
par(mfrow=c(3,2))

plot(matern(varRast[[1]],xyFromCell(varRast,midcell),param=params))
plot(matern(varRast[[1]],xyFromCell(varRast,edgecell),param=params))

plot(varRast[['mid']])
plot(varRast[['edge']])

plot(varRast[['midCor']])
plot(varRast[['edgeCor']])




if(!interactive()){
dev.off()
}

distR = brick(
		distanceFromPoints(varRast,xyFromCell(varRast,midcell)),
		distanceFromPoints(varRast,xyFromCell(varRast,edgecell))
)
names(distR) = c('mid','edge')

if(!interactive()){
pdf("maternGmrfPred.pdf",height=5,width=7)
}
# compare covariance matrix to the matern
xseq = seq(0, 20*params["cellSize"], len=200)
par(mfrow=c(1,1))
plot(xseq, matern(xseq, param=params),
			type = 'l',ylab='cov', xlab='dist',ylim=c(0, params["variance"]),
			main="matern v gmrf")

	# middle cell
points(values(distR[['mid']]), values(varRast[['mid']]), col='red',cex=0.4)
points(values(distR[['mid']]), values(varRast[['midCor']]), col='blue',cex=0.4)


	# edge cells
points(values(distR[['edge']]), values(varRast[['edge']]), col='red',cex=0.6,pch="+")
points(values(distR[['edge']]), values(varRast[['edgeCor']]), col='blue',cex=0.6,pch="+")

	legend("topright", lty=c(1, NA, NA, NA, NA), pch=c(NA, 1, 3, 16, 16),
			col=c('black','black','black','red','blue'),
			legend=c('matern', 'middle','edge','unadj', 'adj')
	)
	if(!interactive()){
		dev.off()	
	}


	covMatMatern = matern(attributes(precMat)$raster, param=params)
	
	prodUncor = covMatMatern %*% precMat
	prodCor = covMatMatern %*% precMatCorr
	
	if(!interactive()){
		pdf("maternGmrfHist.pdf",height=5,width=7)
	}
	par(mfrow=c(2,2))	
	
	hist(Matrix::diag(prodUncor),breaks=40)
	hist(Matrix::diag(prodCor),breaks=40)
	
	hist(prodUncor[lower.tri(prodUncor,diag=FALSE)],breaks=40)	
	hist(prodCor[lower.tri(prodCor,diag=FALSE)],breaks=40)	
	if(!interactive()){
		dev.off()	
	}