library("geostatsp")

# Direct R reference using besselK (old formula)
maternRef = function(x, range = 1, shape = 1, variance = 1) {
	xscale = sqrt(8 * shape) / range
	thisx = abs(x) * xscale
	varscale = variance / (gamma(shape) * 2^(shape - 1))
	res = varscale * thisx^shape * besselK(thisx, shape)
	res[!is.finite(res) & thisx < 1] = variance
	res[!is.finite(res) & thisx >= 1] = 0
	res[x == 0] = variance
	res
}

shapes = c(0.1, 0.5, 0.99, 1, 1.5, 2.5, 4.75, 10)
# distances spanning the short/long branch (logthisx ~ 1.5)
# logthisx = log(d) + 1.5*log(2) + 0.5*log(shape) - log(range)
# with range=1, shape=1: threshold d = exp(1.5 - 1.5*log(2)) ~ 1.59
distances = c(0, 1e-8, 0.01, 0.1, 0.5, 1, 1.5, 2, 5, 10, 50, 100)

maxRelErr = 0
for(shape in shapes) {
	param = c(range = 1, variance = 1.7, shape = shape)
	got = as.numeric(matern(distances, param = param))
	ref = maternRef(distances, range = param["range"],
		shape = param["shape"], variance = param["variance"])
	# relative error where reference is not tiny
	denom = pmax(abs(ref), 1e-300)
	rel = abs(got - ref) / denom
	# absolute tolerance for near-zero values
	absErr = abs(got - ref)
	ok = (rel < 1e-10) | (absErr < 1e-12) | (ref < 1e-300 & got < 1e-300)
	if(!all(ok)) {
		bad = which(!ok)
		stop(sprintf(
			"matern vs besselK mismatch for shape=%g at distances %s\n  got=%s\n  ref=%s\n  rel=%s",
			shape,
			paste(distances[bad], collapse = ","),
			paste(got[bad], collapse = ","),
			paste(ref[bad], collapse = ","),
			paste(rel[bad], collapse = ",")
		))
	}
	maxRelErr = max(maxRelErr, max(rel[ref > 1e-300], 0))
}
cat(sprintf("maternPoint vs besselK: max relative error = %.3e\n", maxRelErr))

# SpatVector path vs matern of Euclidean distances from coordinates
param = c(range = 0.2, shape = 1.5, variance = 1,
	anisoRatio = 1, anisoAngleDegrees = 0)
set.seed(0)
xy = cbind(runif(10), runif(10))
mypoints = vect(xy, crs = "epsg:2000", atts = data.frame(id = 1:10))
myMatern1 = as.matrix(matern(mypoints, param))
dEucl = as.matrix(dist(xy))
myMatern2 = matrix(as.numeric(matern(dEucl, param)), nrow(xy), nrow(xy))
diffMat = abs(myMatern1 - myMatern2)
if(max(diffMat) > 1e-10) {
	stop(sprintf("SpatVector vs Euclidean dist mismatch: max abs err = %g", max(diffMat)))
}
cat(sprintf("SpatVector vs Euclidean dist: max abs error = %.3e\n", max(diffMat)))

# Zero and extreme distances
edge = matern(c(0, 0.001, 100000),
	param = c(range = 1, shape = 1.5, variance = 1))
stopifnot(abs(edge[1] - 1) < 1e-14)
stopifnot(edge[3] < 1e-20)
cat("edge cases OK\n")
