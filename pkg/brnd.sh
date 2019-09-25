R -e "Rcpp::compileAttributes('gpuRandom')"
R -e "devtools::document('gpuRandom')"
R CMD build gpuRandom
tar --list --file=gpuRandom_0.1.tar.gz gpuRandom/inst/doc/
tar --extract --strip=3 --directory=junk --file=gpuRandom_0.1.tar.gz gpuRandom/inst/doc/gpuRn.html
tar --extract --strip=3 --directory=junk --file=gpuRandom_0.1.tar.gz gpuRandom/inst/doc/gpuFisher_sim.html
tar --extract --strip=3 --directory=junk --file=gpuRandom_0.1.tar.gz gpuRandom/inst/doc/matern.html