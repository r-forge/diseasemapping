R -e "Rcpp::compileAttributes('gpuRandom')"
R -e "devtools::document('gpuRandom')"
R CMD build --no-build-vignettes gpuRandom
R CMD INSTALL gpuRandom_0.1.tar.gz
