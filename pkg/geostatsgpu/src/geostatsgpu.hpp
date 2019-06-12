#pragma OPENCL EXTENSION cl_khr_fp64 : enable

//#include <string>
#include <Rcpp.h>
//#include "gpuR/getVCLptr.hpp"

#include "dynVCLMatGeostatsgpu.hpp"
//#include "gpuR/dynVCLMat.hpp"
#include "dynVCLVecGeostatsgpu.hpp"

#include "viennacl/ocl/backend.hpp"
#include "viennacl/linalg/sum.hpp"
#include "viennacl/linalg/lu.hpp"

#define GSL_FLT_EPSILON 1.1920928955078125e-07
#define GSL_DBL_EPSILON 2.2204460492503131e-16
#define GSL_SQRT_DBL_MAX 1.3407807929942596e+154
extern "C" void Rtemme_gamma(double *nu, double * g_1pnu, double * g_1mnu, double *g1, double *g2);

template<typename T> double luT(viennacl::matrix<T> &vclX, viennacl::vector_base<T> &vclD);


template <typename T> T maternClEpsilon();
