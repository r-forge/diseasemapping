
//#include <R_ext/Print.h>
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
//#include "viennacl/linalg/prod.hpp"
//#include "viennacl/linalg/sum.hpp"
//#include "viennacl/linalg/matrix_operations.hpp"


using namespace Rcpp;

extern "C" void Rtemme_gamma(double *nu, double * g_1pnu, double * g_1mnu, double *g1, double *g2);


static const char * kernel_matern =
		"__kernel void matern("
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
		"const double diagVar,"
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
		"const int globalRow = get_global_id(0);"
		"const int Dcol = get_global_id(1);"
		"const int rowHereResult = internalSizeResult*globalRow;"
		"const int rowHereCoords1 = globalRow*internalSizeCoords;"
		"const int rowHereCoords2 = Dcol*internalSizeCoords;"
		"if(globalRow == Dcol){"
		"result[globalRow+rowHereResult] = diagVar;"
		"int Dcol2 = internalSizeCoords*(globalRow+1);"
		//				"result[globalRow+rowHereResult] = globalRow;"
		"result[globalRow+rowHereResult+1] = coords[rowHereCoords1] ;"
		"result[globalRow+rowHereResult+2] = coords[rowHereCoords1 + 1] ;"
		"result[globalRow+rowHereResult+3] = coords[Dcol2] ;"
		"result[globalRow+rowHereResult+4] = coords[Dcol2 + 1] ;"
		"} else if(globalRow < size && Dcol < globalRow){"
		"const int maxIter = 15000;"
		"const double muSq = mu * mu, mup1 = mu+1;"
		"double dist[2], distRotate[2];"
		"double del0, del1, sum0, sum1,fk, pk, qk, hk, ck;"
		"double K_mu, K_mup1, logck;"
		"double K_nu, K_nup1, K_num1, Kp_nu;"
		"int k;"
		"double twoLnHalfX, sigma, half_x_nu, sinhrat;"

		"dist[0] = coords[rowHereCoords1] - coords[rowHereCoords2];"
		"dist[1] = coords[rowHereCoords1 + 1] - coords[rowHereCoords2 + 1];"
		"distRotate[0] = costheta *dist[0] - sintheta * dist[1];"
		"distRotate[1] = sintheta *dist[0] + costheta * dist[1];"
		"const double distSq = distRotate[0]*distRotate[0] + distRotate[1]*distRotate[1]/anisoRatioSq;"
		"const double logthisx = log(distSq- M_LN2) + logxscale;"
		"const double ln_half_x = logthisx - M_LN2;"
		"const double thisx = exp(logthisx);"
		"const double maternBit = varscale + nu * logthisx;"
		"const double expMaternBit = exp(maternBit);"


		// gsl_sf_bessel_K_scaled_temme x < 2
		// gsl_sf_bessel_K_scaled_steed_temme_CF2 x > 2
		"if(logthisx > M_LN2) {"
		"K_nu = 0.0;"
		//	"double bi = 2.0*(1.0 + exp(logthisx));//x);"
		"double logbi = M_LN2 + log1p(thisx);"
		"double bi = exp(logbi);"
		"double di = exp(-logbi);"//1.0/bi;"
		//		  "double delhi = di;" // divide by exp(maternBit)
		"double delhi = exp(-logbi-maternBit);"
		"double hi    = delhi;"// was di, now di/exp(maternBit)

		"double qi   = 0.0;"
		"double qip1 = 1.0;"

		"double ai = -(0.25 - muSq);"
		"double a1 = ai;"
		"double ci = -ai;"
		"double Qi = -ai;"

		"double s = exp(-maternBit) + Qi*delhi;"// was 1 + Qi*delhi
		"  double dels;"
		"  double tmp;"

		"for(k=2; k<=maxIter; k++) {"
		"dels = 0;"
		"tmp = 0;"
		"ai -= 2.0*(k-1);"
		"ci  = -ai*ci/k;"
		"tmp  = (qi - bi*qip1)/ai;"
		"qi   = qip1;"
		"qip1 = tmp;"
		"Qi += ci*qip1;"
		"bi += 2.0;"
		"di  = 1.0/(bi + ai*di);"
		"delhi = (bi*di - 1.0) * delhi;"
		"hi += delhi;"
		"dels = Qi*delhi;"
		"s += dels;"
		"if(fabs(dels/s) < epsilon) break;"
		"}"
		"result[internalSizeResult*Dcol+globalRow] = -k;"

		"hi *= -a1;"

		//		  "*K_nu   = sqrt(M_PI/(2.0*x)) / s;"  sqrt(pi)/2 sqrt(2/x)/s =
		"K_nu= exp(-0.5*ln_half_x)/(M_2_SQRTPI * s);"

		"K_nup1 = K_nu * (mu + thisx + 0.5 - hi*exp(maternBit))/thisx;"
		//		  "Kp_nu  = - *K_nup1 + nu/x * *K_nu;"
		"} else {" // x < 2
		"twoLnHalfX = 2*ln_half_x;"
		"sigma   = - mu * ln_half_x;"
		"half_x_nu = exp(-sigma);"
		"sinhrat = sinh(sigma)/sigma;"

		"fk = sinrat * (cosh(sigma)*g1 - sinhrat*ln_half_x*g2);"
		"pk = 0.5/half_x_nu * g_1pnu;"
		"qk = 0.5*half_x_nu * g_1mnu;"
		"hk = pk;"
		"sum0 = fk*expMaternBit;"
		"sum1 = hk*expMaternBit;"
		"k=0;"
		"logck = maternBit;"
		"del0 = fabs(sum0)+100;"

		"while( (k < maxIter) && ( fabs(del0) > (epsilon * fabs(sum0)) ) ) {"
		"	k++;"
		"   logck += twoLnHalfX - log((double)k);"
		"   ck = exp(logck);"
		"	fk  = (k*fk + pk + qk)/(k*k-muSq);"

		"	del0 = ck * fk;"
		"	sum0 += del0;"

		"	pk /= (k - (mu));"
		"	qk /= (k + (mu));"

		"	hk  = -k*fk + pk;"
		"	del1 = ck * hk;"
		"	sum1 += del1;"
		"}"//while loop
		"result[internalSizeResult*Dcol+globalRow] = k;"
		"}"// end x > 2

		"K_nu   = sum0;"
		"K_nup1 = sum1 * exp( - ln_half_x);"

		"for(k=0; k<nuround; k++) {"
		"	K_num1 = K_nu;"
		"	K_nu   = K_nup1;"
		"	K_nup1 = exp(log(mup1+k) - ln_half_x) * K_nu + K_num1;"
		"}"
		"result[Dcol+rowHereResult] = K_nu;"
		//		"}" // loop col
		"}" // if size
		"}";//function


//[[Rcpp::export]]
void cpp_maternGpu(
		SEXP AR,
		SEXP DR,
		SEXP paramR,
		int max_local_size,
		const int ctx_id) {


	// data
	Rcpp::XPtr<dynVCLMat<double> > pMatA(AR);
	viennacl::matrix_range<viennacl::matrix<double> > A  = pMatA->data();

	Rcpp::XPtr<dynVCLMat<double> > pMatD(DR);
	viennacl::matrix_range<viennacl::matrix<double> > D  = pMatD->data();

	int size=D.size1();
	int internalSizeD = D.internal_size2();
	int internalSizeA = A.internal_size2();

	viennacl::ocl::context ctx(viennacl::ocl::get_context(ctx_id));
	// add kernel to program
	viennacl::ocl::program & my_prog = ctx.add_program(kernel_matern, "kernel_matern");

	// get compiled kernel function
	viennacl::ocl::kernel &matern = my_prog.get_kernel("matern");

	cl_device_id raw_device = ctx.current_device().id();
	cl_kernel raw_kernel = ctx.get_kernel("kernel_matern", "matern").handle().get();

	size_t preferred_work_group_size_multiple;
	cl_int err = clGetKernelWorkGroupInfo(raw_kernel, raw_device,
			CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE,
			sizeof(size_t), &preferred_work_group_size_multiple, NULL);

	max_local_size = roundDown(max_local_size, preferred_work_group_size_multiple);

	Rprintf("s %d %d %d %d %d\n", internalSizeA, internalSizeD, size, max_local_size, preferred_work_group_size_multiple);


	matern.global_work_size(0, internalSizeD);
	matern.global_work_size(1, internalSizeD);

	// set local work sizes
	matern.local_work_size(0, max_local_size);
	matern.local_work_size(1, max_local_size);

	//param is shape range variance nugget anisoRatio anisoAngleRadians anisoAngleDegrees
	double *param = &REAL(paramR)[0];
	int nuround = round(param[0]+0.5);
	double mu = param[0] - nuround;
	double g_1pnu, g_1mnu, g1, g2;
	const double pi_nu = M_PI * mu;
	const double sinrat = (fabs(pi_nu) < GSL_DBL_EPSILON ? 1.0 : pi_nu/sin(pi_nu));
	Rtemme_gamma(&mu, &g_1pnu, &g_1mnu, &g1, &g2);


	// execute kernel
	viennacl::ocl::enqueue(matern(
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
			// nugget + sigmasq
			param[3] + param[2],
			// parameters from bessel temme in gsl
			sinrat, g_1pnu, g_1mnu, g1, g2,
			GSL_DBL_EPSILON /1000,
			A, D));
}

