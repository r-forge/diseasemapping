R -e "Rcpp::compileAttributes('geostatsgpu')"
R -e "devtools::document('geostatsgpu')"
#R CMD build --no-build-vignettes geostatsgpu
#R CMD INSTALL geostatsgpu_0.0.2.tar.gz
R CMD build geostatsgpu
tar --extract --strip=3 --directory=junk --file=geostatsgpu_0.0.2.tar.gz geostatsgpu/inst/doc/matern.html
