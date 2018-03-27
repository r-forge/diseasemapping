
maternKernel = function(
    x, 
    param = c(range = 1, variance = 1, shape = 1), 
    result = vclMatrix(
        data=0, 
        nrow(x), nrow(x),
        type='double',
        ctx_id = x@.context_index),
    type = c("variance", "cholesky", "precision", "inverseCholesky")
) {


	theContstants = .C("maternConstants",
		geostatsp::fillParams(param),
		nuround = -1L,
		mu = -9.9,
		g_1pnu= -9.9, 
		g_1mnu= -9.9,
		g1= -9.9, 
		g2= -9.9,
		sinrat= -9.9,
		sinAngle= -9.9,
		cosAngle= -9.9,
		anisoRatioSq= -9.9,
		varscale= -9.9,
		logxscale= -9.9,
		totalVariance= -9.9) 

gpuSizes = list(size=NA,
	 coordsPadCol = NA, 
	 resultPadRow = NA,
	 resultPadCol = NA)

result = list(
gpuSizes$size,
gpuSizes$coordsPadCol,
gpuSizes$resultPadRow,
gpuSizes$resultPadCol,
theConstants$nu,
theConstants$nuround,
theConstants$mu,
theConstants$costheta,
theConstants$sintheta,
theConstants$anisoRatioSq,
theConstants$varscale,
theConstants$logxscale,
theConstants$diagVar,
theConstants$sinrat,
theConstants$g_1pnu,
theConstants$g_1mnu,
theConstants$g1,
theConstants$g2,
epsilon - 0.000001,
x@coords,
result)



}