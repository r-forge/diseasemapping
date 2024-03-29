# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

backsolveBatchBackend <- function(C, A, B, Cstartend, Astartend, Bstartend, numbatchB, diagIsOne, Nglobal, Nlocal, NlocalCache = 1000L, verbose = 0L) {
    .Call('_gpuRandom_backsolveBatchBackend', PACKAGE = 'gpuRandom', C, A, B, Cstartend, Astartend, Bstartend, numbatchB, diagIsOne, Nglobal, Nlocal, NlocalCache, verbose)
}

cholBatchBackend <- function(A, D, Astartend, Dstartend, numbatchD, Nglobal, Nlocal, NlocalCache) {
    invisible(.Call('_gpuRandom_cholBatchBackend', PACKAGE = 'gpuRandom', A, D, Astartend, Dstartend, numbatchD, Nglobal, Nlocal, NlocalCache))
}

#' Multiply crossproduct matrices
#' 
#' Computes C = t(A) D A
crossprodBatchBackend <- function(C, A, D, invertD, Cstartend, Astartend, Dstartend, Nglobal, Nlocal, NlocalCache, verbose) {
    .Call('_gpuRandom_crossprodBatchBackend', PACKAGE = 'gpuRandom', C, A, D, invertD, Cstartend, Astartend, Dstartend, Nglobal, Nlocal, NlocalCache, verbose)
}

gemmBatchBackend <- function(A, B, C, Arowbatch, Browbatch, Acolbatch, Bcolbatch, need_transpose, Nglobal) {
    .Call('_gpuRandom_gemmBatchBackend', PACKAGE = 'gpuRandom', A, B, C, Arowbatch, Browbatch, Acolbatch, Bcolbatch, need_transpose, Nglobal)
}

#' 
#' Multiplies a rectangular matrix by a rectangular matrix in batches
#'
#' @param C output matrices, stacked row-wise
#' @param A rectangular matrices
#' @param B rectangular matrices 
#' @param Nglobal vector of number of global work items//' @param Nlocal vector of number of local work items//' @param NlocalCache elements in local cache
#' @export
gemmBatch2backend <- function(A, B, C, transposeABC, submatrixA, submatrixB, submatrixC, batches, workgroupSize, NlocalCache, verbose) {
    .Call('_gpuRandom_gemmBatch2backend', PACKAGE = 'gpuRandom', A, B, C, transposeABC, submatrixA, submatrixB, submatrixC, batches, workgroupSize, NlocalCache, verbose)
}

cpp_gpuFisher_test <- function(xR, resultsR, B, streamsR, max_global_size, max_local_size) {
    .Call('_gpuRandom_cpp_gpuFisher_test', PACKAGE = 'gpuRandom', xR, resultsR, B, streamsR, max_global_size, max_local_size)
}

gpuRnBackend <- function(x, streams, max_global_size, random_type) {
    .Call('_gpuRandom_gpuRnBackend', PACKAGE = 'gpuRandom', x, streams, max_global_size, random_type)
}

cpp_gpu_qqnorm <- function(outR, mu, sigma, lowertail, max_global_size, max_local_size) {
    .Call('_gpuRandom_cpp_gpu_qqnorm', PACKAGE = 'gpuRandom', outR, mu, sigma, lowertail, max_global_size, max_local_size)
}

likfitGpu_BackendP <- function(yx, coords, params, boxcox, ssqY, XVYXVX, ssqBetahat, detVar, detReml, jacobian, NparamPerIter, workgroupSize, localSize, NlocalCache, verbose, ssqYX, ssqYXcopy, LinvYX, QinvSsqYx, cholXVXdiag, varMat, cholDiagMat) {
    invisible(.Call('_gpuRandom_likfitGpu_BackendP', PACKAGE = 'gpuRandom', yx, coords, params, boxcox, ssqY, XVYXVX, ssqBetahat, detVar, detReml, jacobian, NparamPerIter, workgroupSize, localSize, NlocalCache, verbose, ssqYX, ssqYXcopy, LinvYX, QinvSsqYx, cholXVXdiag, varMat, cholDiagMat))
}

logfactsumBackend <- function(xR, numWorkItems) {
    .Call('_gpuRandom_logfactsumBackend', PACKAGE = 'gpuRandom', xR, numWorkItems)
}

mat_vec_eledivideBackend <- function(matrixR, rowvectorR, resultR, numWorkItems) {
    invisible(.Call('_gpuRandom_mat_vec_eledivideBackend', PACKAGE = 'gpuRandom', matrixR, rowvectorR, resultR, numWorkItems))
}

matrix_vector_sumBackend <- function(matrixR, rowvectorR, colvectorR, constantR, sumR, numWorkItems) {
    invisible(.Call('_gpuRandom_matrix_vector_sumBackend', PACKAGE = 'gpuRandom', matrixR, rowvectorR, colvectorR, constantR, sumR, numWorkItems))
}

fillParamsExtra <- function(param) {
    invisible(.Call('_gpuRandom_fillParamsExtra', PACKAGE = 'gpuRandom', param))
}

maternBatchBackend <- function(var, coords, param, Nglobal, Nlocal, startrow, numberofrows, verbose = 0L) {
    invisible(.Call('_gpuRandom_maternBatchBackend', PACKAGE = 'gpuRandom', var, coords, param, Nglobal, Nlocal, startrow, numberofrows, verbose))
}

multiplyLowerDiagonalBatchBackend <- function(output, L, D, B, diagIsOne, transformD, Nglobal, Nlocal, NlocalCache) {
    .Call('_gpuRandom_multiplyLowerDiagonalBatchBackend', PACKAGE = 'gpuRandom', output, L, D, B, diagIsOne, transformD, Nglobal, Nlocal, NlocalCache)
}

multiplyDiagonalBatchBackend <- function(C, A, B, inverse, Nglobal, Nlocal) {
    .Call('_gpuRandom_multiplyDiagonalBatchBackend', PACKAGE = 'gpuRandom', C, A, B, inverse, Nglobal, Nlocal)
}

#' Multiply lower triangular matrices
#' 
#' Multiplies a lower triangular matrix by a rectangular matrix
#'
#' @param C output matrices, stacked row-wise
#' @param A lower triangular matrices
#' @param B rectangular matrix or matrices
#' @param Nglobal vector of number of global work items: Drow, Dcol, Dmatrix
#' @param Nlocal vector of number of local work items anything, anything, 1
#' @param NlocalCache elements in local cache
#' @export
multiplyLowerBatchBackend <- function(C, A, B, diagIsOne, Nglobal, Nlocal, NlocalCache) {
    .Call('_gpuRandom_multiplyLowerBatchBackend', PACKAGE = 'gpuRandom', C, A, B, diagIsOne, Nglobal, Nlocal, NlocalCache)
}

createStreamsCpuBackend <- function(n, initial) {
    .Call('_gpuRandom_createStreamsCpuBackend', PACKAGE = 'gpuRandom', n, initial)
}

CreateStreamsGpuBackend <- function(creatorInitialGlobalR, streamsR, keepInitial) {
    invisible(.Call('_gpuRandom_CreateStreamsGpuBackend', PACKAGE = 'gpuRandom', creatorInitialGlobalR, streamsR, keepInitial))
}

