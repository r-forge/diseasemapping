
#include "geostatsgpu.hpp"

using namespace Rcpp;
using namespace viennacl;

// for roundDown
//#include "gpuR/utils.hpp"


double cholGpu(
	viennacl::matrix<double> &x,
	viennacl::vector_base<double> &D,
	viennacl::ocl::kernel &cholKernel
){

	double logdet=0.0; // the result
	unsigned int k;
	int err;

	const unsigned int
		Npad=x.internal_size2(),
		N=x.size2();

	viennacl::ocl::local_mem Dlocal(N*sizeof(cl_double));

	// first column
	viennacl::vector_base<double> firstColX(
			x.handle(), N, 0, 1);
	logdet = firstColX(0);
	D(0) = logdet;
	firstColX *= (1 / logdet);

// remaining columns
	for(k=1;k<N;k++) {
		viennacl::ocl::enqueue(cholKernel(
			x, D, Dlocal, 
			k, N, Npad));
	}

//	logdet = viennacl::linalg::sum(viennacl::linalg::element_log(D));
//	Dlocal = viennacl::linalg::element_log(D);
//	logdet = viennacl::linalg::sum(D);

	return(logdet);
	

}

double cholGpu(
	viennacl::matrix<double> &x,
	viennacl::vector_base<double> &D,
	std::string kernel,
	const int ctx_id,
	int max_global_size
){
	// the context
	viennacl::ocl::context ctx(viennacl::ocl::get_context(ctx_id));
	cl_device_type type_check = ctx.current_device().type();

	// given context but no kernel
	// add kernel to program
	viennacl::ocl::program & my_prog = ctx.add_program(
		kernel,
		"my_kernel");
	// get compiled kernel function
	viennacl::ocl::kernel & cholKernel = my_prog.get_kernel("cholGpu");


	// set global work sizes
    unsigned int M = x.size1();
    unsigned int M_internal = x.internal_size2();

    cl_device_id raw_device = ctx.current_device().id();
    cl_kernel raw_kernel = ctx.get_kernel("my_kernel", "cholGpu").handle().get();    
    size_t preferred_work_group_size_multiple;
        
    cl_int err = clGetKernelWorkGroupInfo(
   		raw_kernel, raw_device,
    	CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, 
        sizeof(size_t), 
        &preferred_work_group_size_multiple, NULL);
        
        if(err != CL_SUCCESS){
            stop("clGetKernelWorkGroupInfo failed");
        }
        
    unsigned int max_local_size = floor(
    	max_global_size/preferred_work_group_size_multiple);

	if(max_local_size <1) max_local_size = 1;

	//	const double workRatio = Ncell/max_local_size;
//	const int workRatioInt = ceil(workRatio);
//	int globalSize = workRatioInt*max_local_size;

	// set work sizes
	// total number of work items
//	cholKernel.global_work_size(0, max_global_size); 
    // number of work items in a group
//	cholKernel.local_work_size(0, max_local_size);

cholKernel.global_work_size(0, 128);//globalSize);
cholKernel.local_work_size(0, 16);

#define DEBUG
#ifdef DEBUG
	Rcout << "global size " << max_global_size << 
		" local size " << max_local_size << 
		" preferred multiple " << 
		preferred_work_group_size_multiple << "\n";
#endif

	double logdet = cholGpu(
		x, D,
		cholKernel);
	return(logdet);

}

//[[Rcpp::export]]
SEXP cpp_cholGpu(
	SEXP xR,
	SEXP DR,
	int max_local_size,
	const int ctx_id,
	SEXP kernelR) {

	double logdet = 0.0;

	// data
	const bool BisVCL=1;
	std::shared_ptr<viennacl::matrix<double> > x = getVCLptr<double>(xR, BisVCL, ctx_id);
	std::shared_ptr<viennacl::vector_base<double> > D = getVCLVecptr<double>(DR, BisVCL, ctx_id);

    std::string kernel = as<std::string>(kernelR);


	logdet = cholGpu(
		*x, *D, kernel,
		ctx_id, 
		max_local_size);

	return(Rcpp::wrap(logdet));	

}