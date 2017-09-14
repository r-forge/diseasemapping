
//#include <R_ext/Print.h>
#include <Rmath.h>

#define GSL_DBL_EPSILON 2.2204460492503131e-16
#define GSL_SQRT_DBL_MAX 1.3407807929942596e+154

// make link to /usr/lib/x86_64-linux-gnu/libOpenCL.so

//#include <RcppEigen.h>
#include <Rcpp.h>

#include "viennacl/ocl/backend.hpp"

#include "gpuR/utils.hpp"
#include "gpuR/getVCLptr.hpp"


extern "C" void Rtemme_gamma(double *nu, double * g_1pnu, double * g_1mnu, double *g1, double *g2);

using namespace Rcpp;

//[[Rcpp::export]]
void cpp_maternGpu(
		SEXP AR,
		SEXP DR,
		SEXP paramR,
		SEXP sourceCode_,
//		SEXP cholSourceCode_,
//		int type,
		int max_local_size,
		const int ctx_id) {

	std::string my_kernel = as<std::string>(sourceCode_);
	viennacl::ocl::context ctx(viennacl::ocl::get_context(ctx_id));

	// data
	const bool BisVCL=1;
	viennacl::matrix<double> *vclA, *vclD;
	vclA = getVCLptr<double>(AR, BisVCL, ctx_id);
	vclD = getVCLptr<double>(DR, BisVCL, ctx_id);

	int sizeD1=vclD->size1(),sizeD2=vclD->size2(),sizeA1=vclA->size1(),sizeA2=vclA->size2();
	int iSizeD1=vclD->internal_size1(),iSizeD2=vclD->internal_size2();
	int iSizeA1=vclA->internal_size1(),iSizeA2=vclA->internal_size2();


	// add kernel to program
	viennacl::ocl::program & my_prog = ctx.add_program(my_kernel, "my_kernel");

	// get compiled kernel function
	viennacl::ocl::kernel & maternKernel = my_prog.get_kernel("maternCL");

	cl_device_type type_check = ctx.current_device().type();


	/*
	 * n by n, n*(n-1)/2 total
	 * Dindex = ( n*(n-1) - Dcol*(Dcol-1) )/2 + Drow
	 * Dcol*Dcol - Dcol - 2*Drow + 2*Dindex - n*(n-1) = 0
	 * Dcol = (1 + sqrt(1 + 4*n*n-4*n - 8*Dindex - 8*Drow)) )/2
	 */
	int Ntriangle = sizeD1 * (sizeD1 - 1)/2;

	if(max_local_size > sizeD1) max_local_size = sizeD1;

	// set global work sizes
	const double workRatio = sizeD1/max_local_size;
	const int workRatioInt = ceil(workRatio);
	int globalSize = workRatioInt*max_local_size;
	maternKernel.global_work_size(0, globalSize);
	maternKernel.global_work_size(1, globalSize);

	// set local work sizes
	maternKernel.local_work_size(0, max_local_size);
	maternKernel.local_work_size(1, max_local_size);

//	Rprintf("a1 %d a2 %d d1 %d d2 %d a1 %d a2 %d d1 %d d2 %d s %d m %d r %d\n", sizeA1, sizeA2, sizeD1, sizeD2,iSizeA1, iSizeA2, iSizeD1, iSizeD2, max_local_size, Ntriangle,  roundDown(Ntriangle, sizeD1));


	//param is shape range variance nugget anisoRatio anisoAngleRadians anisoAngleDegrees
	double *param = &REAL(paramR)[0];
	int nuround = round(param[0]+0.5);
	double mu = param[0] - nuround;
	double g_1pnu, g_1mnu, g1, g2;
	const double pi_nu = M_PI * mu;
	const double sinrat = (fabs(pi_nu) < GSL_DBL_EPSILON ? 1.0 : pi_nu/sin(pi_nu));
	Rtemme_gamma(&mu, &g_1pnu, &g_1mnu, &g1, &g2);

	// execute kernel

	viennacl::ocl::enqueue(maternKernel(
			sizeA1, iSizeA2, iSizeD1, iSizeD2,
			// nuround mu cos theta, sin theta
			param[0], nuround, mu, cos(param[5]), sin(param[5]),
			// parameters from matern.c in geostatsp
			// anisoRatioSq
			(param[4])*(param[4]),
			// varscale
			log(param[2])  - Rf_lgammafn(param[0]) -  (param[0]-1)*M_LN2,
			// logxscale
			1.5 * M_LN2 +   0.5 * log(param[0]) - log(param[1]),
			// nugget + sigmasq
			param[3] + param[2],
			// parameters from bessel temme in gsl
			sinrat, g_1pnu, g_1mnu, g1, g2,
			GSL_DBL_EPSILON /1000,
			*vclA, *vclD));


}

