#' @title LU decomposition
#'
#' @description LU decomposition using ViennaCL
#' 
#' @param x matrix
#' @examples
#' library('geostatsgpu')
#' gpuR::setContext(grep('gpu', listContexts()$device_type) [1])
#' v1 = matrix(1,10,10)
#' diag(v1) = 2
#' v2 <- gpuR::vclMatrix(v1, type='float')
#' v3 = gpuR::deepcopy(v2)
#' system.time(res1 <- luGpu(v2))
#' system.time(res2 <- chol(v1))
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
