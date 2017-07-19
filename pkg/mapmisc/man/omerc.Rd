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
    preserve=NULL) 
    ocea(x, angle, flip=FALSE)
    moll(x=0, angle=NULL, flip=FALSE)
}

\arguments{
	\item{x}{A SpatialP* object or a vector of length 2 giving the 
	centroid of the projection. }
\item{angle}{angle of rotation or vector of angles}
\item{post}{post-projection angle rotation}
\item{flip}{post-projection flipping of coordinates}
\item{preserve}{A \code{SpatialPoints} object, the resulting projection is scaled
so as to preserve the distances between these points as best as possible. }
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
	An object of class \code{\link[sp]{CRS}}.
}	

\references{
\url{http://en.wikipedia.org/wiki/Space-oblique_Mercator_projection}
}
\seealso{\code{\link[sp]{CRS}},\code{\link[sp]{spTransform}}}

\examples{

omerc(c(10,50), angle=c(0,45,80))

data('netherlands')

if(require('rgdal', quietly=TRUE)){
nldUtm = spTransform(nldCities,
	omerc(nldCities, angle=0))
projection(nldUtm)

map.new(nldUtm)
text(nldUtm,labels=nldUtm$name)
scaleBar(nldUtm, 'topright')

nldRot = spTransform(nldCities,
omerc(nldCities, angle=seq(25,45,by=5))
)
projection(nldRot)

map.new(nldRot)
text(nldRot,labels=nldRot$name)
scaleBar(nldRot, 'topright')
insetMap(nldRot, 'bottomright',map=world)

}
\dontrun{

	nldOceaCrs = ocea(nldCities)
	nldOcea = spTransform(nldCities, nldOceaCrs)
	map.new(nldOcea, buffer=2000)
	text(nldOcea,labels=nldOcea$name)
	scaleBar(nldOcea, 'topright')

	map.new(nldCities)
	plot(nldTiles,add=TRUE)
	text(nldCities,labels=nldCities$name)

	tilesRot = openmap(nldRot)
	map.new(nldRot)
	plot(tilesRot,add=TRUE)
	text(nldRot,labels=nldRot$name)
	
	tilesUtm = openmap(nldUtm)
	map.new(nldUtm)
	plot(tilesUtm,add=TRUE)
	text(nldUtm,labels=nldUtm$name)

}

}