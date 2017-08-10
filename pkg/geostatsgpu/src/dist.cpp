
#include <R_ext/Print.h>
#include <Rmath.h>

#define GSL_DBL_EPSILON 2.2204460492503131e-16
#define GSL_SQRT_DBL_MAX 1.3407807929942596e+154

// make link to /usr/lib/x86_64-linux-gnu/libOpenCL.so

// Use OpenCL with ViennaCL
#define VIENNACL_WITH_OPENCL 1

#include "gpuR/dynVCLMat.hpp"
#include "Rcpp.h"
#include <RcppEigen.h>
#include "gpuR/utils.hpp"


// ViennaCL headers
#include "viennacl/ocl/device.hpp"
#include "viennacl/ocl/platform.hpp"
#include "viennacl/matrix.hpp"
#include "viennacl/linalg/prod.hpp"
#include "viennacl/linalg/sum.hpp"
//#include "viennacl/linalg/matrix_operations.hpp"


using namespace Rcpp;

RcppExport void Rtemme_gamma(double *nu, double * g_1pnu, double * g_1mnu, double *g1, double *g2);

//http://viennacl.sourceforge.net/doc/custom-kernels_8cpp_source.html
static const char * my_compute_program_bessel =
		"__kernel void besselk(\n"
		"          __global const float * dist,\n"
		"          __global const float * nu, \n"
		"          __global float * result,\n"
		"          unsigned int size) \n"
		"{ \n"
		"		int Dn, e10;"
		"		int k;"
		"		int stat_iter;"
		//		"		float mu = *muR;"
		"		float K_mu, K_mup1, Kp_mu;"
		"		float K_nu, K_nup1, K_num1, Kp_nu;"
		"		float pi_nu = M_PI * (*nu);"
		"		float half_x, ln_half_x, half_x_nu, sigma, sinhrat, ex;"
		"		float sum0, sum1;"
		"		float fk, pk, qk, hk, ck;"
		// maxIter, sinrat, pi_nu, mu, nuround
		//		"		double g_1pnu, g_1mnu, g1, g2;"
		"		double del0;"
		"		double del1;"

		"  for (unsigned int i = get_global_id(0); i < size; i += get_global_size(0))\n"
		"k=0;"
		"ex = exp(dist[i]);"
		"half_x    = 0.5 * x[D];"
		"ln_half_x = log(half_x);"
		"half_x_nu = exp((*nu)*ln_half_x);"
		"sigma   = - (*nu) * ln_half_x;"

		"    result[i] = dist[i] * (*nu);\n"
		"};\n";

static const char * kernel_matern =
		"__kernel void anisoDist("
		"const int size,"
		"const int internalSizeCoords,"
		"const int internalSizeResult,"
		"const double nu,"
		"const int nuround,"
		"const double mu,"
		"const double costheta,"
		"const double sintheta,"
		"const double anisoRatioSq,"
		"const double varscale,"
		"const double logxscale,"
		"const double sinrat,"
		"const double g_1pnu,"
		"const double g_1mnu,"
		"const double g1,"
		"const double g2,"
		"const double epsilon,"
		"__global const double *coords,"
		"__global double *result)"
		"{"
		// Get the index of the elements to be processed
		"const int maxIter = 15000;"
		"const int globalRow = get_global_id(0);"
		"const int globalCol = get_global_id(1);"
		"const int rowHereResult = internalSizeResult*globalRow;"
		"const int rowHereCoords1 = internalSizeCoords*globalRow;"
		"const int rowHereCoords2 = internalSizeCoords*globalCol;"

		"if((globalRow < size) && (globalCol < globalRow)){"

		"double dist[2], distRotate[2], distSq;"
		"const double muSq = mu * mu, mup1 = mu+1;"
		"double del0, del1, sum0, sum1,fk, pk, qk, hk, ck;"
		"double K_mu, K_mup1, sumk;"
		"double K_nu, K_nup1, K_num1, Kp_nu;"
		"int k;"

		"dist[0] = coords[rowHereCoords1] - coords[rowHereCoords2];"
		"dist[1] = coords[rowHereCoords1 + 1] - coords[rowHereCoords2 + 1];"
		"distRotate[0] = costheta *dist[0] - sintheta * dist[1];"
		"distRotate[1] = sintheta *dist[0] + costheta * dist[1];"
		"distSq = distRotate[0]*distRotate[0] + distRotate[1]*distRotate[1]/anisoRatioSq;"

		"const double logthisx = 0.5 * log(distSq) + logxscale;"
		"const double ln_half_x = logthisx - M_LN2;"
		"const double twoLnHalfX = 2*ln_half_x;"
		"const double maternBit = varscale + nu * logthisx;"
		"const double sigma   = - mu * ln_half_x;"
		"const double half_x_nu = exp(mu*ln_half_x);"//exp(-sigma);"
		"const double sinhrat = sinh(sigma)/sigma;"
		"const double half_x = exp(ln_half_x);" // get rid?

		"fk = sinrat * (cosh(sigma)*g1 - sinhrat*ln_half_x*g2);"
		"pk = 0.5/half_x_nu * g_1pnu;"
		"qk = 0.5*half_x_nu * g_1mnu;"
		"ck = exp(maternBit);"
		"hk = pk;"
		"sum0 = fk*exp(maternBit);"
		"sum1 = hk*exp(maternBit);"
		"k=0;"
		"sumk = 0;"//maternBit;"
		"del0 = fabs(sum0);"

		"while( (k < maxIter) && ( fabs(del0) > (epsilon * fabs(sum0)) ) ) {"
		"	k++;"
		"   sumk += log((double)k);"
		"	fk  = (k*fk + pk + qk)/(k*k-muSq);"
		"	ck *= half_x*half_x/k;"
		// log(ck) += 2*log(x) - log(4) - log(k)"
		// log(ck) = k*log4 - sumk +2*k*log(x)
//		"   ck = exp(maternBit + 2* k * logthisx - sumk);"
		"	del0 = ck * fk;"
		"	sum0 += del0;"

		"	pk /= (k - (mu));"
		"	qk /= (k + (mu));"

		"	hk  = -k*fk + pk;"
		"	del1 = ck * hk;" // = ck * (pk - k fk)
		"	sum1 += del1;"
		"}"//while loop

		"K_nu   = sum0;"// * exp(maternBit);"
		"K_nup1 = sum1 * exp( - ln_half_x);"

		"for(k=0; k<nuround; k++) {"
		"	K_num1 = K_nu;"
		"	K_nu   = K_nup1;"
			// does this need modifying if we're doing the matern?
		"	K_nup1 = exp(log(mup1+k) - ln_half_x) * K_nu + K_num1;"
		"}"
		"result[globalCol*internalSizeResult+globalRow] = k;"

		"result[globalCol+rowHereResult] = K_nu;"
		"}" // if size
		"}";//function

RcppExport SEXP cpp_gpuMatrix_custom_P(
				SEXP AR,
				SEXP DR,
				SEXP paramR,
				SEXP max_local_sizeR,
				SEXP ctx_idR) {

	//    Rcpp::traits::input_parameter< int >::type max_local_size(max_local_sizeR);
	int max_local_size = as<int>(max_local_sizeR);
	Rcpp::traits::input_parameter< const int >::type ctx_id(ctx_idR);

	// data
	Rcpp::XPtr<dynVCLMat<double> > pMatA(AR);
	viennacl::matrix_range<viennacl::matrix<double> > A  = pMatA->data();

	Rcpp::XPtr<dynVCLMat<double> > pMatD(DR);
	viennacl::matrix_range<viennacl::matrix<double> > D  = pMatD->data();

	int size=D.size1();
	int internalSizeD = D.internal_size1();
	int internalSizeA = A.internal_size1();

	viennacl::ocl::context ctx(viennacl::ocl::get_context(ctx_id));
	// add kernel to program
	viennacl::ocl::program & my_prog = ctx.add_program(kernel_matern, "kernel_matern");

	// get compiled kernel function
	viennacl::ocl::kernel & anisoDist = my_prog.get_kernel("anisoDist");

	cl_device_id raw_device = ctx.current_device().id();
	cl_kernel raw_kernel = ctx.get_kernel("kernel_matern", "anisoDist").handle().get();

	size_t preferred_work_group_size_multiple;
	cl_int err = clGetKernelWorkGroupInfo(raw_kernel, raw_device,
			CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE,
			sizeof(size_t), &preferred_work_group_size_multiple, NULL);

	max_local_size = floor(max_local_size / preferred_work_group_size_multiple);
	anisoDist.global_work_size(0, internalSizeD);
	anisoDist.global_work_size(1, internalSizeD);

	// set local work sizes
	anisoDist.local_work_size(0, max_local_size);
	anisoDist.local_work_size(1, max_local_size);


	//param is shape range variance nugget anisoRatio anisoAngleRadians anisoAngleDegrees
	double *param = &REAL(paramR)[0];
	int nuround = round(param[0]+0.5);
	double mu = param[0] - nuround;
	double g_1pnu, g_1mnu, g1, g2;
	const double pi_nu = M_PI * mu;
	const double sinrat = (fabs(pi_nu) < GSL_DBL_EPSILON ? 1.0 : pi_nu/sin(pi_nu));
	Rtemme_gamma(&mu, &g_1pnu, &g_1mnu, &g1, &g2);

//	Rprintf("s %f %d %f %f %f %f %f \n", mu, nuround, sinrat, g_1pnu, g_1mnu, g1, g2);

	// execute kernel
	viennacl::ocl::enqueue(anisoDist(
			size, internalSizeA, internalSizeD,
			// nuround mu cos theta, sin theta
			param[0], nuround, mu, cos(param[5]), sin(param[5]),
			// parameters from matern.c in geostatsp
			// anisoRatioSq
			(param[4])*(param[4]),
			// varscale
			log(param[2])  - Rf_lgammafn(param[0]) -  (param[0]-1)*M_LN2,
			// logxscale
			1.5 * M_LN2 +   0.5 * log(param[0]) - log(param[1]),
			// parameters from bessel temme in gsl
			sinrat, g_1pnu, g_1mnu, g1, g2,
			GSL_DBL_EPSILON * 0.5,
			A, D));
	return R_NilValue;
}

