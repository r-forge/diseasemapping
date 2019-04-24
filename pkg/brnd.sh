R -e "Rcpp::compileAttributes('clRNG')"
R -e "devtools::document('clRNG')"
R CMD build clRNG
tar --extract --strip=3 --directory=junk --file=clRNG_0.0.1.tar.gz clRNG/inst/doc/random.html
