#pragma OPENCL EXTENSION cl_khr_fp64 : enable

#include <clRNG/clRNG.h>
#include <Rcpp.h>
#include <string>

#include <dynMatrix/dynVCLMatGeostatsgpu.hpp>
#include <dynMatrix/dynVCLVecGeostatsgpu.hpp>
#include "viennacl/linalg/sum.hpp"
#include "viennacl/ocl/backend.hpp"

//from typeDef.cpp
template <typename T> std::string openclTypeString();
template <typename T> int sizeOfReal();






#include <clRNG/mrg31k3p.h>

// clRNG -> Matrix
void convertclRngMat(clrngMrg31k3pStream* streams, Rcpp::IntegerMatrix result);
//matrix ->clRNG streams
void convertMatclRng(Rcpp::IntegerMatrix Sin, clrngMrg31k3pStream* streams);

template <typename T> std::string mrg31k3pTypeString();









//////////////////////////////////////////////////////////////////////////////////////////////////////
//#include "viennacl/linalg/lu.hpp"

extern "C" void Rtemme_gamma(double *nu, double * g_1pnu, double * g_1mnu, double *g1, double *g2);

//template<typename T> double luT(viennacl::matrix<T> &vclX, viennacl::vector_base<T> &vclD);

template <typename T> T maternClEpsilon();


















