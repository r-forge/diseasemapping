#' Convert rasters to lgcp objects
#' 
#' @description Convert a RasterLayer to spatialAtRisk
#' @param X an object of class \code{\link[raster]{RasterLayer-class}}
#' @param ... not used
#' 
#' @return A \code{fromXYZ} object produced by 
#' \code{\link[lgcp]{spatialAtRisk}} 
#' 
#' @examples
#' library('raster')
#' x = raster(matrix(1:12, 4,3,byrow=TRUE))/100
#' sum(values(x))*prod(res(x))
#' par(mfrow=c(1,2))
#' plot(x)
#' \dontrun{
#' library('lgcp')
#' x2 = spatialAtRisk(x)
#' sum(x2$Zm)*prod(res(x))
#' plot(x2)
#' }
#' @export
spatialAtRisk.RasterLayer = function(X, ...) {
	
	# replace NA's with zero
	X = raster::subs(X, data.frame(NA,0), subsWithNA=FALSE)
	
	# convert to xyz list
	Xlist = list(
			X=sort(raster::xFromCol(X)), 
			Y = sort(raster::yFromRow(X)), 
			Z=as.matrix(raster::t(raster::flip(X, direction='y'))))
	# convert to fromXYZ object
	lgcp::spatialAtRisk(Xlist)
	
}
