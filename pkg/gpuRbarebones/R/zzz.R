#' @useDynLib gpuRbarebones
#' @importFrom Rcpp evalCpp
#' @importFrom utils strOptions packageVersion tail

.onLoad <- function(libname, pkgname) {
    options(gpuR.print.warning=TRUE)
    options(gpuR.default.type = "double")
    # options(gpuR.default.device.type = "gpu")
}

.onAttach <- function(libname, pkgname) {
 
    if (!identical(Sys.getenv("APPVEYOR"), "True")) {
        # initialize contexts
    #      initContexts()

        packageStartupMessage(paste0("gpuRbarebones ", packageVersion('gpuRbarebones')))
    }
}

.onUnload <- function(libpath) {
    options(gpuR.print.warning=NULL)
    options(gpuR.default.type = NULL)
    # options(gpuR.default.device.type = NULL)
}
