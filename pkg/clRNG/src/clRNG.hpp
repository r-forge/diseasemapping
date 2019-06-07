#pragma OPENCL EXTENSION cl_khr_fp64 : enable

#include <clRNG/clRNG.h>
#include <Rcpp.h>


#include "dynVCLMatGeostatsgpu.hpp"
#include "dynVCLVecGeostatsgpu.hpp"
#include "viennacl/linalg/sum.hpp"
#include "viennacl/ocl/backend.hpp"

