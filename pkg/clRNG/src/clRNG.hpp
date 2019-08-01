//#pragma OPENCL EXTENSION cl_khr_fp64 : enable

#include <clRNG/clRNG.h>
#include <Rcpp.h>
#include <string>


#include "dynVCLMatGeostatsgpu.hpp"
#include "dynVCLVecGeostatsgpu.hpp"
#include "viennacl/linalg/sum.hpp"
#include "viennacl/ocl/backend.hpp"


template <typename T> std::string openclTypeString();
template <typename T> int sizeOfReal();
