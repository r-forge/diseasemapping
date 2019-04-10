R -e "Rcpp::compileAttributes('gpuRbarebones')"
R -e "pkgbuild::compile_dll('gpuRbarebones')"
R -e "devtools::document('gpuRbarebones')"
 R CMD build gpuRbarebones
# R CMD INSTALL geostatsgpu_0.0.2.tar.gz
#R CMD build gpuRbarebones
tar --extract --strip=3 --directory=junk --file=gpuRbarebones_0.0.1.tar.gz gpuRbarebones/inst/doc/gpuRbarebones.pdf
