// Copyright (c) 2016 Wladimir J. van der Laan
// Distributed under the MIT software license.
// Based on an example from the OpenCL cookbook.


//#include "gpuRbarebones/windows_check.hpp"
//#include "gpuRbarebones/utils.hpp"
#include <CL/cl.h>
//#include "gpuR/cl_helpers.hpp"
// Use OpenCL with ViennaCL
#define VIENNACL_WITH_OPENCL 1

#include <Rcpp.h>

using namespace Rcpp;


#include "viennacl/ocl/device.hpp"
#include "viennacl/ocl/platform.hpp"
#include "viennacl/ocl/backend.hpp"


struct platform_data_item {
    int id;
    char* name;
};

struct platform_data_item platform_data_items[] = {
    { CL_PLATFORM_PROFILE, "Profile"},
    { CL_PLATFORM_VERSION, "Version"},
    { CL_PLATFORM_NAME,    "Name"},
    { CL_PLATFORM_VENDOR,  "Vendor"},
    { CL_PLATFORM_EXTENSIONS, "Extensions"},
};

#define ARRAYLEN(array)     (sizeof(array)/sizeof((array)[0]))

//' @title Current Device Information
//' @description Check current device information
//' @return list containing
//' @return \item{device}{Character string of device name}
//' @return \item{device_index}{Integer identifying device}
//' @return \item{device_type}{Character string identifying device type (e.g. gpu)}
//' @export
// [[Rcpp::export]]
SEXP currentDevice()
{
    std::string device_type;

    

    cl_device_type check = viennacl::ocl::current_device().type(); 

    if(check & CL_DEVICE_TYPE_CPU){
    device_type = "cpu";
    }else if(check & CL_DEVICE_TYPE_GPU){
    device_type = "gpu";
    }else if(check & CL_DEVICE_TYPE_ACCELERATOR){
    device_type = "accelerator";
    }else{
    Rcpp::Rcout << "device found: " << std::endl;
    Rcpp::Rcout << check << std::endl;
    throw Rcpp::exception("unrecognized device detected");

    }
    
    return List::create(Named("device") = wrap(viennacl::ocl::current_context().current_device().name()),
                        Named("device_type") = wrap(device_type));

}


//' @title Detect Number of Platforms
//' @description Find out how many OpenCL enabled platforms are available.
//' @return An integer value representing the number of platforms available.
//' @details see https://laanwj.github.io/2016/05/06/opencl-ubuntu1604.html
//' @seealso \link{detectGPUs}
//' @export
// [[Rcpp::export]]
int detectPlatforms2() {
 
    int i, j;
    char value[1024];
    char data[1024];
    size_t retsize;
    size_t valueSize;
    cl_uint platformCount = 0;
    cl_platform_id platforms[128];
    cl_uint deviceCount;
    cl_device_id devices[128];
    cl_uint maxComputeUnits;
    // get all platforms
//    clGetPlatformIDs(0, platforms, &platformCount);

    if (clGetPlatformIDs(0, NULL, &platformCount) != CL_SUCCESS) {
        stop("Unable to get platform IDs\n");
    }
//    platforms = (cl_platform_id*) malloc(sizeof(cl_platform_id) * platformCount);
    if (clGetPlatformIDs(platformCount, platforms, NULL) != CL_SUCCESS) {
        stop("Unable to get platform IDs\n");
    }
 
    for (i = 0; i < platformCount; i++) {
        Rcout << i+1 << ". Platform\n";

        for (int j=0; j<ARRAYLEN(platform_data_items); ++j) {
            if (clGetPlatformInfo(platforms[i], platform_data_items[j].id, sizeof(data), data, &retsize) != CL_SUCCESS || retsize == sizeof(data)) {
                Rcpp::Rcerr << "Unable to get platform" << platform_data_items[j].name << "\n";
            }
            Rcout << platform_data_items[j].name << ": " << data << "\n";
        }
 
        // get all devices
        if (clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, 0, NULL, &deviceCount) != CL_SUCCESS) {
            Rcpp::Rcerr << "Unable to get device IDs for platform " << platforms[i] << "\n";
        }
//        devices = (cl_device_id*) malloc(sizeof(cl_device_id) * deviceCount);
        if (clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, deviceCount, devices, NULL) != CL_SUCCESS) {
            Rcpp::Rcerr << "Unable to get device IDs for platform " << platforms[i] << "\n";
        }

 
        // for each device print critical attributes
        for (j = 0; j < deviceCount; j++) {
 
            // print device name
            clGetDeviceInfo(devices[j], CL_DEVICE_NAME, 0, NULL, &valueSize);
            clGetDeviceInfo(devices[j], CL_DEVICE_NAME, valueSize, value, NULL);
            Rcout << " " << j+1 << ". Device: " << value << "\n";

 
            // print hardware device version
            clGetDeviceInfo(devices[j], CL_DEVICE_VERSION, 0, NULL, &valueSize);

            clGetDeviceInfo(devices[j], CL_DEVICE_VERSION, valueSize, value, NULL);
            Rcout << " " << j+1 << "." << 1 << " Hardware version: " <<  value<< "\n";

 
            // print software driver version
            clGetDeviceInfo(devices[j], CL_DRIVER_VERSION, 0, NULL, &valueSize);

            clGetDeviceInfo(devices[j], CL_DRIVER_VERSION, valueSize, value, NULL);
            Rcout << " " << j+1 << "." << 2 << " Software version: "  << value << "\n";

 
            // print c version supported by compiler for device
            clGetDeviceInfo(devices[j], CL_DEVICE_OPENCL_C_VERSION, 0, NULL, &valueSize);

            clGetDeviceInfo(devices[j], CL_DEVICE_OPENCL_C_VERSION, valueSize, value, NULL);
            Rcout << " " <<j+1 << "."<< 3 << " OpenCL C version: " << value << "\n";
 
            // print parallel compute units
            clGetDeviceInfo(devices[j], CL_DEVICE_MAX_COMPUTE_UNITS,
                    sizeof(maxComputeUnits), &maxComputeUnits, NULL);
            Rcout << " " << j+1 << "."<< 4 << " Parallel compute units: " << maxComputeUnits << "\n";
 
        }
 

 
    }


    return 0;
 
}