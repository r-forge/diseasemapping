.onLoad <- function(libname, pkgname) {
  # Set OMP_NUM_THREADS only when running CMD check, examples, tests, or vignettes
  is_check <- isTRUE(as.logical(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))) ||
              isTRUE(as.logical(Sys.getenv("R_PACKAGE_CHECKING"))) ||
              any(grepl("CheckExEnv", search())) ||
              !is.null(Sys.getenv("R_CMD_CHECK", unset = NULL))

  if (is_check) {
    Sys.setenv(OMP_NUM_THREADS = "2")
  }
}