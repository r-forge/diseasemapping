
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

//' @title Detect Number of Platforms
//' @description Find out how many OpenCL enabled platforms are available.
//' @return An integer value representing the number of platforms available.
//' @seealso \link{detectGPUs}
//' @export
// [[Rcpp::export]]
SEXP detectPlatforms()
{
    platforms_type platforms = viennacl::ocl::get_platforms();
    
    return wrap(platforms.size());
}

//// [[Rcpp::export]]
//void setPlatform(SEXP platform_idx_)
//{
//    unsigned int platform_idx = as<unsigned int>(platform_idx_);
//    typedef std::vector< viennacl::ocl::platform > platforms_type;
//    
//    // get platforms
//    platforms_type platforms = viennacl::ocl::get_platforms();
//    unsigned int platforms_size = platforms.size();
//    
//    if(platform_idx > platforms_size){
//        stop("Platform index out of bounds");
//    }
//    
//    // subtract one for zero indexing    
//    // set platform
//    long id = 0;
//    viennacl::ocl::set_context_platform_index(id, platform_idx - 1);
//}

//' @title Return Current Platform
//' @description Find out which platform is currently in use
//' @return \item{platform}{Name of the current platform}
//' @return \item{platform_index}{Index of current platform}
//' @seealso \link{detectPlatforms}
//' @export
// [[Rcpp::export]]
SEXP currentPlatform()
{
    // get current platform index
    int plat_idx = viennacl::ocl::current_context().platform_index();
    
    // get platforms
    platforms_type platforms = viennacl::ocl::get_platforms();
    
//    return wrap(plat_idx + 1);
//    return wrap(platforms[plat_idx].info());
    
    return List::create(Named("platform") = wrap(platforms[plat_idx].info()),
                        Named("platform_index") = wrap(plat_idx + 1));
}

//List platformNames()
//{
//    
//    std::cout << platforms[0].info();
//}

//std::vector<std::string> split(const std::string &s, char delim) {
//    std::vector<std::string> elems;
//    split(s, delim, elems);
//    return elems;
//}


// [[Rcpp::export]]
List cpp_platformInfo(SEXP platform_idx_)
{
    cl_int err;
    
    // get platforms
    platforms_type platforms = viennacl::ocl::get_platforms();
    
    // subtract one for zero indexing
    unsigned int platform_idx = as<unsigned int>(platform_idx_) - 1;
    
    viennacl::ocl::platform vcl_platform = platforms[platform_idx];
    
    cl_platform_id platform_id = vcl_platform.id();
    
    char platformName[1024];
    char platformVendor[1024];
    char platformVersion[1024];
    char platformExtensions[2048];
        
    err = clGetPlatformInfo(platform_id, 
                            CL_PLATFORM_NAME,
                            sizeof(platformName),
                            platformName,
                            NULL
    );
    
    if(err != CL_SUCCESS){
        Rcpp::stop("Acquiring platform ID failed");
    }
    
    err = clGetPlatformInfo(platform_id, 
                            CL_PLATFORM_VENDOR,
                            sizeof(platformName),
                            platformVendor,
                            NULL
    );
    
    if(err != CL_SUCCESS){
        Rcpp::stop("Acquiring platform vendor failed");
    }
    
    err = clGetPlatformInfo(platform_id, 
                            CL_PLATFORM_VERSION,
                            sizeof(platformName),
                            platformVersion,
                            NULL
    );
    
    if(err != CL_SUCCESS){
        Rcpp::stop("Acquiring platform version failed");
    }
    
    err = clGetPlatformInfo(platform_id, 
                            CL_PLATFORM_EXTENSIONS,
                            sizeof(platformName),
                            platformExtensions,
                            NULL
    );
    
    if(err != CL_SUCCESS){
        Rcpp::stop("Acquiring platform extensions failed");
    }
    
    // Convert char arrays to string
    std::string platformNameStr(platformName);
    std::string platformVendorStr(platformVendor);
    std::string platformVersionStr(platformVersion);

    // Split extensions to a vector
    std::vector<std::string> extensionsVector;
    std::string platformExtensionsStr(platformExtensions);
    extensionsVector = split(platformExtensionsStr, ' ');

    return List::create(Named("platformName") = platformNameStr,
                        Named("platformVendor") = platformVendorStr,
                        Named("platformVersion") = platformVersionStr,
                        Named("platformExtensions") = extensionsVector
                        );
}



//' @title Details of platforms
//' @description Details of OpenCL enabled platforms which are available.
//' @return A list
//' @seealso \link{detectGPUs}
//' @export
// [[Rcpp::export]]
List platformsAvailable()
{
    cl_int err;
    char value[1024];
    size_t valueSize, retsize;

    // get platforms
    cl_uint deviceCount, platformCount;
    cl_platform_id platforms[128];
    cl_device_id devices[128];

    int Dplatform, Ddevice, Nplatforms, Ndevices;

    int Nitems = 7;


    if (clGetPlatformIDs(0, NULL, &platformCount) != CL_SUCCESS) {
        stop("Unable to get platform IDs\n");
    }
    if (clGetPlatformIDs(platformCount, platforms, NULL) != CL_SUCCESS) {
        stop("Unable to get platform IDs\n");
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
