
#include "gpuRbarebones/windows_check.hpp"
#include "gpuRbarebones/utils.hpp"

// Use OpenCL with ViennaCL
#define VIENNACL_WITH_OPENCL 1

// ViennaCL headers
#include "viennacl/ocl/backend.hpp"
#include "viennacl/ocl/platform.hpp"
#include "viennacl/detail/matrix_def.hpp"

#include <Rcpp.h>

using namespace Rcpp;

typedef std::vector< viennacl::ocl::platform > platforms_type;



//' @title Details of platforms
//' @description Details of OpenCL enabled platforms which are available.
//' @return A list
//' @seealso \link{detectGPUs}
//' @export
// [[Rcpp::export]]
List platformsAvailable() {

    char value[1024];
    size_t valueSize, retsize;

    // get platforms
    cl_uint deviceCount, platformCount;
    cl_platform_id platforms[128];
    cl_device_id devices[128];

    int Dplatform, Ddevice, Nplatforms, Ndevices;

    int Nitems = 7;


    if (clGetPlatformIDs(0, NULL, &platformCount) != CL_SUCCESS) {
        warning("Unable to get platform IDs\n");
        platformCount = 0;
    }
    if (clGetPlatformIDs(platformCount, platforms, NULL) != CL_SUCCESS) {
        warning("Unable to get platform IDs (1)\n");
        platformCount = 0;
    }

    Nplatforms = platformCount;
    CharacterVector platformVector(5);
    platformVector.names() = CharacterVector::create("NAME", "VENDOR", "VERSION", "PROFILE", "EXTENSIONS");

    List result(Nplatforms);

    for(Dplatform=0; Dplatform < Nplatforms; Dplatform++) {
        if (clGetDeviceIDs(platforms[Dplatform], CL_DEVICE_TYPE_ALL, 0, NULL, &deviceCount) != CL_SUCCESS) {
            Rcpp::Rcerr << "Unable to get device IDs for platform " << platforms[Dplatform] << "\n";
        }
        if (clGetDeviceIDs(platforms[Dplatform], CL_DEVICE_TYPE_ALL, deviceCount, devices, NULL) != CL_SUCCESS) {
            Rcpp::Rcerr << "Unable to get device IDs for platform " << platforms[Dplatform] << "\n";
        }

        clGetPlatformInfo(platforms[Dplatform], CL_PLATFORM_NAME, sizeof(value), value, &retsize);
        platformVector[0] = value;
        clGetPlatformInfo(platforms[Dplatform], CL_PLATFORM_VENDOR, sizeof(value), value, &retsize);
        platformVector[1] = value;
        clGetPlatformInfo(platforms[Dplatform], CL_PLATFORM_VERSION,  sizeof(value), value, &retsize);
        platformVector[2] = value;
        clGetPlatformInfo(platforms[Dplatform], CL_PLATFORM_PROFILE,  sizeof(value), value, &retsize);
        platformVector[3] = value;
//        clGetPlatformInfo(platforms[Dplatform], CL_PLATFORM_EXTENSIONS,  sizeof(value), value, &retsize);
//        platformVector[4] = value;

        Ndevices = deviceCount;

        CharacterMatrix deviceMat(Nitems, Ndevices);
        rownames(deviceMat) = CharacterVector::create(
            "NAME", "VENDOR", "VERSION", "PROFILE", 
            "OPENCL_C_VERSION", "DRIVER_VERSION", "EXTENSIONS");


        for(Ddevice = 0; Ddevice < Ndevices; Ddevice++) {
            clGetDeviceInfo(devices[Ddevice], CL_DEVICE_VENDOR, 0, NULL, &valueSize);
            clGetDeviceInfo(devices[Ddevice], CL_DEVICE_VENDOR, valueSize, value, NULL);
            deviceMat[1 + Nitems * Ddevice] = value;

            clGetDeviceInfo(devices[Ddevice], CL_DEVICE_NAME, 0, NULL, &valueSize);
            clGetDeviceInfo(devices[Ddevice], CL_DEVICE_NAME, valueSize, value, NULL);
            deviceMat[0+ Nitems * Ddevice] = value;
            clGetDeviceInfo(devices[Ddevice], CL_DEVICE_VERSION, 0, NULL, &valueSize);
            clGetDeviceInfo(devices[Ddevice], CL_DEVICE_VERSION, valueSize, value, NULL);
            deviceMat[2+ Nitems * Ddevice] = value;
            clGetDeviceInfo(devices[Ddevice], CL_DEVICE_PROFILE, 0, NULL, &valueSize);
            clGetDeviceInfo(devices[Ddevice], CL_DEVICE_PROFILE, valueSize, value, NULL);
            deviceMat[3+ Nitems * Ddevice] = value;
            clGetDeviceInfo(devices[Ddevice], CL_DEVICE_OPENCL_C_VERSION, 0, NULL, &valueSize);
            clGetDeviceInfo(devices[Ddevice], CL_DEVICE_OPENCL_C_VERSION, valueSize, value, NULL);
            deviceMat[4+ Nitems * Ddevice] = value;
            clGetDeviceInfo(devices[Ddevice], CL_DRIVER_VERSION, 0, NULL, &valueSize);
           clGetDeviceInfo(devices[Ddevice], CL_DRIVER_VERSION, valueSize, value, NULL);
            deviceMat[5+ Nitems * Ddevice] = value;
//            clGetDeviceInfo(devices[Ddevice], CL_DEVICE_EXTENSIONS, valueSize, value, NULL);
//            clGetDeviceInfo(devices[Ddevice], CL_DEVICE_EXTENSIONS, 0, NULL, &valueSize);
//            deviceMat[6+ Nitems * Ddevice] = value;

        }


        result[Dplatform] = List::create(
            Named("platform") = Rcpp::clone(platformVector),
            Named("devices") = Rcpp::clone(deviceMat)
            );
    }

    return result;
}
