#include <clRNG/mrg31k3p.h>
#include "clRNG.hpp"
#include <string>
#include "viennacl/ocl/backend.hpp"


// clRNG -> Matrix
void convertclRngMat(clrngMrg31k3pStream* streams, Rcpp::IntegerMatrix result);
//matrix ->clRNG streams
void convertMatclRng(Rcpp::IntegerMatrix Sin, clrngMrg31k3pStream* streams);