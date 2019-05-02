R -e "Rcpp::compileAttributes('gpuMatrix')"
R -e "devtools::document('gpuMatrix')"
R CMD build gpuMatrix
tar --extract --strip=3 --directory=junk --file=gpuMatrix_0.0.1.tar.gz gpuMatrix/inst/doc/lowerMult.html
