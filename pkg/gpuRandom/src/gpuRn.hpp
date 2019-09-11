//#pragma OPENCL EXTENSION cl_khr_fp64 : enable

#include <Rcpp.h>
#include <string>


#include "viennacl/linalg/sum.hpp"
#include "viennacl/ocl/backend.hpp"

#include "dynVCL/dynVCLMatGeostatsgpu.hpp"
#include "dynVCL/dynVCLVecGeostatsgpu.hpp"


#include <clRNG/mrg31k3p.h>
#include <clRNG/clRNG.h>




// clRNG -> Matrix
void convertclRngMat(clrngMrg31k3pStream* streams, Rcpp::IntegerMatrix result);
//matrix ->clRNG streams
void convertMatclRng(Rcpp::IntegerMatrix Sin, clrngMrg31k3pStream* streams);


// in gpuRn.cpp
template <typename T> std::string mrg31k3pTypeString();

//from typeDef.cpp
template <typename T> std::string openclTypeString();
template <typename T> int sizeOfReal();
template <typename T> T maternClEpsilon();