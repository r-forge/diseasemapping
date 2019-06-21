#include <clRNG/mrg31k3p.h>
#include "clRNG.hpp"


// clRNG -> Matrix
void convertclRngMat(clrngMrg31k3pStream* streams, Rcpp::IntegerMatrix result);
//matrix ->clRNG streams
void convertMatclRng(Rcpp::IntegerMatrix Sin, clrngMrg31k3pStream* streams);

template <typename T> std::string mrg31k3pTypeString();