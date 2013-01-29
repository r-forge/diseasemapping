\name{maternGmrfPrec}
\alias{maternGmrfPrec}
\title{Precision matrix for a Matern spatial correlation}
\description{
Produces the precision matrix for a Gaussian random field on a regular square lattice, using a Markov random field approximation. 
}

\usage{ 
maternGmrfPrec(Nx, Ny=Nx, range, scale=1/range, var=1, prec=1/var, kappa=2,cellSize=1) 
}

\arguments{
\item{Nx}{Number of grid cells in the x direction}
\item{Ny}{Grid cells in the y direction, defaults to \code{Nx} for a square grid}
\item{range}{Range parameter for the correlation function, in units of distance rather than as grid cells.  Ignored if \code{scale} is specified}
\item{scale}{Scale parameter for the correlation function.}
\item{var}{Variance parameter, ignored if \code{predc} is specified}
\item{prec}{Precision parameter}
\item{kappa}{Roughness parameter of the Matern correlation}
\item{cellSize}{The size of the grid cells.}
}
\value{
A precision matrix, contained in a sparse matrix object. }
\details{ The numbering of cells
is consistent with INLA and the \code{raster} package.  Cell 1 is the top left cell, with cell 2 being the cell to the right and numbering
continuing row-wise.
}

\examples{
	range = 3
	cellSize=2
	#Nx = Ny = 60
	Nx=20;Ny=25
	kappa=1
	var=2^2
	precMat =maternGmrfPrec(Nx, Ny, range=range, var=var, 
		kappa=kappa,cellSize=cellSize) 
	
	mycell = Nx*(floor(Ny/2)-1) + floor(Nx/2) # the middle cell
	a=matrix(precMat[,mycell], Ny, Nx, byrow=TRUE)
 	a[seq(floor(Ny/2)-3, by=1, len=7), seq(floor(Nx/2)-3, by=1, len=7)]
	
	varMat = solve(precMat)
	b=matrix(varMat[,mycell], Ny, Nx, byrow=TRUE)
	
	xseq = seq(-(Nx-1)*(cellSize/2) + cellSize/2, by=cellSize, len=Nx)
	plot(xseq, b[floor(Ny/2),], ylim=c(0, max(c(var, max(b)))))
	xseq = seq(min(xseq), max(xseq), len=1000)
	xseq2 = abs(xseq/range) 
	lines(xseq, var*(2^(1-kappa)/gamma(kappa)) * xseq2^kappa * besselK(xseq2, nu=kappa), col='red')
	
\dontrun{	
	cseq = seq(1,len=n,by=cellSize)
	coords = c(outer(1i*cseq, cseq, FUN="+"))
	distmat = outer(coords, coords, FUN="-") 
	distmat = abs(distmat)/range
	covmat = (2^(1-kappa)/gamma(kappa)) * distmat^kappa * besselK(distmat, nu=kappa)
	diag(covmat) = 1
	covmat = covmat * var
	
	exactPrecMat = solve(covmat)

	c=matrix(exactPrecMat[,mycell], Ny, Nx, byrow=TRUE)[aseq, aseq]
	round(c*100)/100

	prodmat =  covmat %*% precMat
	round(100*as.matrix(prodmat[seq(floor(Ny/2)-3, by=1, len=7), seq(floor(Nx/2)-3, by=1, len=7)]))/100	
}

}
