#' @title LU decomposition
#'
#' @description LU decomposition using ViennaCL
#' 
#' @param x matrix
#' @examples
#' library('geostatsgpu')
#' v1 <- maternGpu(
#' 	as.matrix(expand.grid(seq(1,by=1,len=23), seq(201,by=1,len=14))),
#' 	c(shape=4.5, range=1.5, variance = 2, nugget = 0, anisoRatio = 2, anisoAngleDegrees = 12.5))
#' v2 = as.matrix(v1) 
#' v3 = deepcopy(v1)
#' system.time(res1 <- luGpu(v1))
#' system.time(res2 <- chol(v2))
#' system.time(res3 <- cholGpu(v3))
#' 
#' @export
luGpu = function(x, D) {

	if(missing(D)) {
		D = vclVector(0, nrow(x),
			type=typeof(x),
			ctx_id = x@.context_index)
	}

	logDet=cpp_lu(
		x,
		D)

	list(L=x, D=D, logDet = logDet)

}
