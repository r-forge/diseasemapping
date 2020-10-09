maternGpuParam = function(x, type='double') {
	x = geostatsp::fillParam(x)
	paramsBatch = vclMatrix(cbind(x, 
		matrix(0, nrow(x), 22-ncol(x))),type=type)
	paramsBatch
}