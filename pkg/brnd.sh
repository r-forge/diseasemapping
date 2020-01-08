rm junk/fisher_sim.html
rm gpuRandom_*.gz
R -e "Rcpp::compileAttributes('gpuRandom')"
R -e "devtools::document('gpuRandom')"
R CMD build gpuRandom
tar --list --file=gpuRandom_0.1.tar.gz gpuRandom/inst/doc/
tar --extract --strip=3 --directory=junk --file=gpuRandom_0.1.tar.gz gpuRandom/inst/doc/random_number.html
tar --extract --strip=3 --directory=junk --file=gpuRandom_0.1.tar.gz gpuRandom/inst/doc/fisher_sim.html
tar --extract --strip=3 --directory=junk --file=gpuRandom_0.1.tar.gz gpuRandom/inst/doc/matern.html
tar --extract --strip=3 --directory=junk --file=gpuRandom_0.1.tar.gz gpuRandom/inst/doc/qqnorm.html