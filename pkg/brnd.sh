R -e "Rcpp::compileAttributes('gpuRandom')"
R -e "devtools::document('gpuRandom')"
R CMD build --no-build-vignettes gpuRandom
R CMD INSTALL gpuRandom_0.1.tar.gz




# R -e "Rcpp::compileAttributes('geostatsp')"
# R -e "devtools::document('geostatsp')"
# R CMD build --no-build-vignettes geostatsp
# R CMD INSTALL geostatsp_1.9.0.tar.gz