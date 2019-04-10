#' @useDynLib gpuRbarebones, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom utils strOptions packageVersion tail



.onAttach <- function(libname, pkgname) {

        packageStartupMessage(paste0("gpuRbarebones ", packageVersion('gpuRbarebones')))

}