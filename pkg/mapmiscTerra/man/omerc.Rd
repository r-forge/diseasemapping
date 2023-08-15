\name{omerc}
\alias{omerc}
\alias{moll}
\alias{ocea}
\title{Oblique Mercator, Cylindrical, and Mollweide projections}
\description{
Defines an appropriate Oblique Mercator, Oblique Cylindrical Equal Area, and Mollweide projections
for a supplied Spatial object
}
\usage{
	omerc(x, angle,  
	post=c('none', 'north', 'wide','tall'),
    preserve=NULL, ellipse=TRUE) 
    ocea(x, angle, flip=FALSE)
    moll(x=0, angle=NULL, flip=FALSE)
}

\arguments{
	\item{x}{A \code{\link[terra]{SpatVector}} object or a vector of length 2 giving the 
	centroid of the projection. }
\item{angle}{angle of rotation or vector of angles}
\item{post}{post-projection angle rotation}
\item{flip}{post-projection flipping of coordinates}
\item{preserve}{A \code{\link[terra]{SpatVector}} object, the resulting projection is scaled
so as to preserve the distances between these points as best as possible. }
\item{ellipse}{compute projection region and areas to crop when projecting.}
}

\details{
With \code{omerc}, an Oblique Mercator map projection is produced which warps the world onto a cylinder, with the north-south axis 
rotated by the specified angle.  If \code{angle} is a vector, the optimal 
angle for reducing the size
of the bounding box is returned.   

If \code{post = 'north'}, an inverse rotation will preserve the north direction at the origin.

If \code{post = 'wide'}, an inverse rotation
makes the smallest possible bounding box which is wider than tall.  

If \code{post = 'tall'}, the bounding box is taller than it is wide

If \code{post} is numeric, it specifies an angle for inverse rotation.

\code{ocea} produces an Oblique Cylindrical Equal Area projection and \code{moll} a Mollweide projections
}


\value{
	An object of class \code{\link[terra]{crs}}.
}	

\references{
\url{https://en.wikipedia.org/w/index.php?title=Space-oblique_Mercator_projection}
}

\seealso{
  \code{\link[terra]{crs}},
  \code{\link[terra]{project}}
}

\examples{

data('worldMap')
worldMap = terra::unwrap(worldMap)


myProj = omerc(c(-100,-70), angle=-45)
crs(myProj, proj=TRUE)

plot(project(worldMap, crsLL))
plot(attributes(myProj)$crop, col='red', add=TRUE)



}