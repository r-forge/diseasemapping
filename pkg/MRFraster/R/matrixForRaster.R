matrixForRaster <-
function(x) {	
	t(x)[ncol(x):1,]
}
