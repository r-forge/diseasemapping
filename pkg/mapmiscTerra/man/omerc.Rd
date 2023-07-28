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

data('netherlands')
nldCities = terra::unwrap(nldCities)
nldTiles = terra::unwrap(nldTiles)


nldUtm = project(nldCities, omerc(nldCities, angle=0))
cat(crs(nldUtm, proj=TRUE))
nldTiles2 = project(nldTiles, crs(nldUtm), method='near')


map.new(nldTiles2)
plot(nldTiles2, add=TRUE)
points(nldUtm)
text(nldUtm,labels=nldUtm$name, pos=4)
scaleBar(nldUtm, 'topright')
insetMap(crs=nldUtm, pos='bottomright', cropInset = 20*1000*c(2,2,1,3), zoom=7, 
	inset=c(0.01, 0.01), width=0.25)
insetMap(nldUtm, 'right', cropInset = 500*1000, zoom=4, width=0.25, inset = 0)


	nldRot = project(nldCities, omerc(nldCities, angle=seq(25,45,by=5)))
	cat(crs(nldRot, proj=TRUE))


# testing the wrapping of polygons with an ugly projection
crsMollRot = moll(nldCities, angle = 70)
map.new(crsMollRot, col='lightblue')
plot(wrapPoly(worldMap, crsMollRot), col='red', add=TRUE)

plot(project(worldMap, crsLL))
plot(attributes(crsMollRot)$crop, add=TRUE, col='red')

crsMollRot = moll(worldMap[grep("Sierra Leone", worldMap$NAME)] , angle = -45)
map.new(crsMollRot, col='lightblue')
plot(wrapPoly(worldMap, crsMollRot), col='red', add=TRUE)

plot(project(worldMap, crsLL))
plot(attributes(crsMollRot)$crop, add=TRUE, col='red')

\dontrun{

map.new(nldRot)
text(nldRot,labels=nldRot$name)
scaleBar(nldRot, 'topright')
insetMap(nldRot, 'bottomright')

crsMoll = moll()
worldMap2 = project(worldMap, crsMoll)
plot(attributes(crsMoll)$ellipse, col='lightblue')
plot(worldMap2, add=TRUE, col='yellow')

crsMoll3 = moll(-100)
worldMap3 = wrapPoly(worldMap, crsMoll3)
plot(attributes(crsMoll3)$ellipse, col='lightblue')
plot(worldMap3, add=TRUE, col='yellow')


crsO = ocea(c(0,1), flip='wsu')
worldMapO = wrapPoly(worldMap, crsO)
plot(attributes(crsO)$ellipse, col='lightblue')
plot(worldMapO, add=TRUE, col='yellow')
}


}