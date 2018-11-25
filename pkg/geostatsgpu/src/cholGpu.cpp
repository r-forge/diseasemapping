
#include "geostatsgpu.hpp"

using namespace Rcpp;
using namespace viennacl;

double cholGpu(
	viennacl::matrix<double> &x,
	viennacl::vector_base<double> &D,
	viennacl::ocl::kernel &cholKernel
){

	double logdet=0.0; // the result
	unsigned int k;

	const unsigned int 
		iSize1=x.internal_size1(),
		size1=x.size1();


	viennacl::vector<double> Dlocal(iSize1);

	for(k=0;k<size1;k++) {
	viennacl::ocl::enqueue(cholKernel(
		x, D, Dlocal, k, iSize1, size1))
	}

//	logdet = viennacl::linalg::sum(element_log(D));
	logdet = viennacl::linalg::sum(D);

	return(logdet);
	

}

double cholGpu(
	viennacl::matrix<double> &x,
	viennacl::vector_base<double> &DofLDL,
	char kernel,
	const int ctx_id,
	int max_local_size,
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
	viennacl::ocl::kernel & cholKernel = my_prog.get_kernel("maternCL");

	// set global work sizes
	const unsigned int 
		size1=x.size1(),
		size2=x.size2(),
		Ncell = size1 * (size1 - 1)/2;

	if(max_local_size > Ncell) max_local_size = Ncell;

	const double workRatio = Ncell/max_local_size;
	const int workRatioInt = ceil(workRatio);
	int globalSize = workRatioInt*max_local_size;

	// set work sizes
	cholKernel.global_work_size(0, globalSize);
	cholKernel.local_work_size(0, max_local_size);

	double logdet = maternGpuVcl(
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
	const char kernel) {

	double logdet = 0.0;

	// data
	const bool BisVCL=1;
	std::shared_ptr<viennacl::matrix<double> > x = getVCLptr<double>(xR, BisVCL, ctx_id);
	std::shared_ptr<viennacl::vector_base<double> > D = getVCLVecptr<double>(DR, BisVCL, ctx_id);

	logdet = cholGpu(
		*x, *D, kernel,
		ctx_id, 
		max_local_size);

	return(Rcpp::wrap(logdet));	

}