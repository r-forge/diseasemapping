#pragma OPENCL EXTENSION cl_khr_fp64 : enable

//#include <string>
#include <Rcpp.h>

#include "dynVCLMatGpuMatrix.hpp"
#include "dynVCLVecGpuMatrix.hpp"

#include "viennacl/ocl/backend.hpp"
#include "viennacl/linalg/sum.hpp"


