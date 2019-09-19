//#pragma OPENCL EXTENSION cl_khr_fp64 : enable

#include <clRNG/clRNG.h>
#include <Rcpp.h>
#include <string>


#include "vclMatrix/dynVCLMatGeostatsgpu.hpp"
#include "vclMatrix/dynVCLVecGeostatsgpu.hpp"
#include "viennacl/linalg/sum.hpp"
#include "viennacl/ocl/backend.hpp"

//from typeDef.cpp
template <typename T> std::string openclTypeString();
template <typename T> int sizeOfReal();


// clRNG -> Matrix
void convertclRngMat(clrngMrg31k3pStream* streams, Rcpp::IntegerMatrix result);
//matrix ->clRNG streams
void convertMatclRng(Rcpp::IntegerMatrix Sin, clrngMrg31k3pStream* streams);

template <typename T> std::string mrg31k3pTypeString();