//#pragma OPENCL EXTENSION cl_khr_fp64 : enable

#include <clRNG/mrg31k3p.h>
#include <clRNG/clRNG.h>

#include <Rcpp.h>
#include <string>


#include "dynVCLMatGeostatsgpu.hpp"
#include "dynVCLVecGeostatsgpu.hpp"
#include "viennacl/linalg/sum.hpp"
#include "viennacl/ocl/backend.hpp"

#include <CL/mrg31k3pkernelStringSeparate.hpp>


// clRNG -> Matrix
void convertclRngMat(clrngMrg31k3pStream* streams, Rcpp::IntegerMatrix result);
//matrix ->clRNG streams
void convertMatclRng(Rcpp::IntegerMatrix Sin, clrngMrg31k3pStream* streams);

template <typename T> std::string mrg31k3pTypeString();

//from typeDef.cpp
template <typename T> std::string openclTypeString();
template <typename T> int sizeOfReal();







