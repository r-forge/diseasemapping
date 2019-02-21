#include "geostatsgpu.hpp"
#include "viennacl/ocl/backend.hpp"

using namespace Rcpp;
using namespace viennacl;

//[[Rcpp::export]] 
Rcpp::List gpuNlocal(
	std::string kernel,
	std::string functionName,
	int ctx_id
	){

//    std::string kernel = as<std::string>(kernelR_);
//    std::string functionName = as<std::string>(functionNameR_);


	// the context
	viennacl::context ctx(viennacl::ocl::get_context(ctx_id));
	viennacl::ocl::context ctx2(viennacl::ocl::get_context(ctx_id));
    unsigned int raw_device = ctx.opencl_context().current_device_id();
//	cl_device_type type_check = ctx.current_device().type();
	viennacl::ocl::device working_device = ctx.opencl_context().devices()[raw_device];

    int
    	localMem = working_device.local_mem_size(),
    	sizeOfDouble = sizeof(cl_double),
    	maxWorkgroupSize = working_device.max_work_group_size();

	// add kernel to program

	viennacl::ocl::program & my_prog = ctx2.add_program(
		kernel, "my_kernel");

	// get compiled kernel function
	viennacl::ocl::kernel & cholKernel = 
		my_prog.get_kernel(
			functionName
			);

    cl_kernel raw_kernel = ctx2.get_kernel(
    	"my_kernel", 
    	functionName
    	).handle().get();

    size_t preferred_work_group_size_multiple;
    cl_device_id thedeviceid;    

    unsigned int err = clGetKernelWorkGroupInfo(
   		raw_kernel, thedeviceid,
    	CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, 
        sizeof(size_t), 
        &preferred_work_group_size_multiple, NULL);
        
        if(err != CL_SUCCESS){
            stop("clGetKernelWorkGroupInfo failed");
        }

	return Rcpp::List::create(
		Rcpp::Named("localWorkgroupSize") = Rcpp::wrap(preferred_work_group_size_multiple),
		Rcpp::Named("localMemory") = Rcpp::wrap(localMem),
		Rcpp::Named("sizeOfDouble") = Rcpp::wrap(sizeOfDouble),
		Rcpp::Named("maxWorkgroupSize") = Rcpp::wrap(maxWorkgroupSize)
		);
}
