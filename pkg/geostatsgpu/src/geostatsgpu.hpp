#pragma OPENCL EXTENSION cl_khr_fp64 : enable

//#include <string>
#include <Rcpp.h>

#include "dynVCLMatGeostatsgpu.hpp"
#include "dynVCLVecGeostatsgpu.hpp"

#include "viennacl/ocl/backend.hpp"
#include "viennacl/linalg/sum.hpp"
#include "viennacl/linalg/lu.hpp"


extern "C" void Rtemme_gamma(double *nu, double * g_1pnu, double * g_1mnu, double *g1, double *g2);

template<typename T> double luT(viennacl::matrix<T> &vclX, viennacl::vector_base<T> &vclD);


template <typename T> int sizeOfReal();
template <typename T> std::string openclTypeString();
template <typename T> T maternClEpsilon();