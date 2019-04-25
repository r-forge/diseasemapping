#pragma OPENCL EXTENSION cl_khr_fp64 : enable

//#include <string>
#include <Rcpp.h>

#include "dynVCLMatGpuMatrix.hpp"
#include "dynVCLVecGpuMatrix.hpp"

#include "viennacl/ocl/backend.hpp"
#include "viennacl/linalg/sum.hpp"


template <typename T> 
std::string openclTypeString() {
  return("undefined");}

template <> std::string openclTypeString<double>(){
  std::string result = "double";
  return(result);
}

template <> std::string openclTypeString<float>(){
  std::string result = "float";
  return(result);
}

template <> std::string openclTypeString<int>(){
  std::string result = "int";
  return(result);
}