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
