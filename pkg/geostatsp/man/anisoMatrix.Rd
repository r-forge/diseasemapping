\name{anisoMatrix}
\alias{anisoMatrix}
Geometric Anisotropy Matrix
\title{
Calculates a 2-D matrix for geometric anisotropy
}
\usage{
anisoMatrix(angle, scale)
}
\arguments{
  \item{angle}{
	angle of rotaiton, in radians
}
  \item{scale}{
	scale for x and y axes after rotation
}
}
\details{
Produces matrix A such that
	(xnew, ynew) = A (x,y)	
}
\value{
A 2 by 2 matrix
}
\author{
Patrick Brown
}


\seealso{
\code{\link{GaussRF}}
}
\examples{

myAnisoMat = anisoMatrix(-pi/8,c(2, 16))

mymodel=list("$", var=2,   
		A=myAnisoMat,
		list("whittle", 2))
		
myraster = raster(nrows=41,ncols=41,xmn=0,xmx=10,ymn=0,ymx=10)

myGRF=GaussRF(myraster, model=mymodel)

plot(myGRF)


}