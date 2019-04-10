#include "gpuRbarebones/windows_check.hpp"
// #include "gpuR/cl_helpers.hpp"

// Use OpenCL with ViennaCL
#define VIENNACL_WITH_OPENCL 1
//#define VIENNACL_DEBUG_ALL 1

// ViennaCL headers
#include "viennacl/ocl/device.hpp"
#include "viennacl/ocl/platform.hpp"
#include "viennacl/ocl/backend.hpp"

#include <Rcpp.h>

// using namespace cl;
using namespace Rcpp;

//' @title Initiate a single context
//' @description Initiate a viennaCL context on the supplied platform and device
//' @return An integer value representing the context id
//' @export
// [[Rcpp::export]]
SEXP initSingleContext(int platformId, int deviceId, int contextId){

    unsigned int plat_idx = platformId-1, gpu_idx = deviceId-1, id = contextId-1;

    // get platforms
    typedef std::vector< viennacl::ocl::platform > platforms_type;
    platforms_type platforms = viennacl::ocl::get_platforms();

    std::vector< viennacl::ocl::device > devices;
    devices = platforms[plat_idx].devices(CL_DEVICE_TYPE_ALL);


    viennacl::ocl::set_context_platform_index(id, plat_idx);
    viennacl::ocl::setup_context(id, devices[gpu_idx]);

    
    viennacl::ocl::switch_context(id);
    
    return Rcpp::wrap((int) (id+1));
}




//' @title Info for an OpenCL Context
//' @description Provide information on a specified context
//' @return data frame with one row and elements
//' @return \item{context}{Integer identifying context}
//' @return \item{platform}{Character string listing OpenCL platform}
//' @return \item{platform_index}{Integer identifying platform}
//' @return \item{device}{Character string listing device name}
//' @return \item{device_index}{Integer identifying device}
//' @return \item{device_type}{Character string labeling device (e.g. gpu)}
//' @export
// [[Rcpp::export]]
DataFrame contextInfo(int contextId) {

    unsigned int id = contextId-1;
    int num_contexts=1;


    Rcpp::IntegerVector   context_index(num_contexts);
    Rcpp::IntegerVector   platform_index(num_contexts);
    Rcpp::CharacterVector device_name(num_contexts);
    

    viennacl::ocl::switch_context(id);
    context_index[0] = id + 1;

    viennacl::vcl_size_t plat_idx = viennacl::ocl::current_context().platform_index();
    platform_index[0] = (int) (plat_idx)+1;

    device_name[0] = viennacl::ocl::current_context().current_device().name();



    
    return Rcpp::List::create(
        Rcpp::Named("context") = context_index,
//                  Rcpp::Named("platform") = platform_name,
                  Rcpp::Named("platform_index") = platform_index,
                  Rcpp::Named("device") = device_name);
}

