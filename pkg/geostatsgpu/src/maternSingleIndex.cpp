
//#include <R_ext/Print.h>
#include <Rmath.h>

#define GSL_DBL_EPSILON 2.2204460492503131e-16
#define GSL_SQRT_DBL_MAX 1.3407807929942596e+154

// make link to /usr/lib/x86_64-linux-gnu/libOpenCL.so

//#include <RcppEigen.h>
#include <Rcpp.h>
#include <RcppEigen.h>

#include "viennacl/linalg/lu.hpp"
#include "viennacl/ocl/backend.hpp"

#include "gpuR/utils.hpp"
#include "gpuR/getVCLptr.hpp"


extern "C" void Rtemme_gamma(double *nu, double * g_1pnu, double * g_1mnu, double *g1, double *g2);

using namespace Rcpp;

//[[Rcpp::export]]
void cpp_maternGpuSingleIndex(
	SEXP varR,
	SEXP DofLDLR,
		SEXP XYR, // solve for Lt b = XY
		SEXP crossprodR, //bt b
		SEXP coordsR,
		SEXP paramR,
		const int type, // 2 cholesky 3 inversecholesky, 4 inverse, 5 solve for b
		const int upper,
		SEXP sourceCode_,
		int max_local_size,
		const int ctx_id) {

	std::string my_kernel = as<std::string>(sourceCode_);
	viennacl::ocl::context ctx(viennacl::ocl::get_context(ctx_id));

	// add kernel to program
	viennacl::ocl::program & my_prog = ctx.add_program(my_kernel, "my_kernel");

	// get compiled kernel function
	viennacl::ocl::kernel & maternKernel = my_prog.get_kernel("maternSingleIndexCL");

	cl_device_type type_check = ctx.current_device().type();

	// data
	const bool BisVCL=1;
	unsigned int D;

	std::shared_ptr<viennacl::matrix<double> > vclVar = getVCLptr<double>(varR, BisVCL, ctx_id);
	std::shared_ptr<viennacl::matrix<double> > vclCoords = getVCLptr<double>(coordsR, BisVCL, ctx_id);
	std::shared_ptr<viennacl::matrix<double> > vclDofLDL = getVCLptr<double>(DofLDLR, BisVCL, ctx_id);


	const unsigned int 
	sizeCoords1=vclCoords->size1(),
	sizeCoords2=vclCoords->size2(),
	sizeVar1=vclVar->size1(),
	sizeVar2=vclVar->size2(),
	iSizeCoords1=vclCoords->internal_size1(),
	iSizeCoords2=vclCoords->internal_size2(),
	iSizeVar1=vclVar->internal_size1(),
	iSizeVar2=vclVar->internal_size2(),
	Ncell = sizeVar1 * (sizeVar1 - 1)/2;

	if(max_local_size > Ncell) max_local_size = Ncell;

	// set global work sizes
	const double workRatio = Ncell/max_local_size;
	const int workRatioInt = ceil(workRatio);
	int globalSize = workRatioInt*max_local_size;


	// set work sizes
	maternKernel.global_work_size(0, globalSize);
	maternKernel.local_work_size(0, max_local_size);

//	Rprintf("a1 %d a2 %d d1 %d d2 %d a1 %d a2 %d d1 %d d2 %d s %d m %d r %d\n", sizeA1, sizeA2, sizeD1, sizeD2,iSizeA1, iSizeA2, iSizeD1, iSizeD2, max_local_size, Ntriangle,  roundDown(Ntriangle, sizeD1));


	//param is shape range variance nugget anisoRatio anisoAngleRadians anisoAngleDegrees
	double *param = &REAL(paramR)[0];
	int nuround = round(param[0]+0.5);
	const unsigned int maxIter = 1500;
	double mu = param[0] - nuround;
	double g_1pnu, g_1mnu, g1, g2;
	const double muSq = mu*mu, 
	mup1 = mu + 1.0,
	pi_nu = M_PI * mu;
	const double sinrat = (fabs(pi_nu) < GSL_DBL_EPSILON ? 1.0 : pi_nu/sin(pi_nu));

	Rtemme_gamma(&mu, &g_1pnu, &g_1mnu, &g1, &g2);

	// execute kernel

	viennacl::ocl::enqueue(maternKernel(
		Ncell, iSizeCoords2, iSizeVar1, iSizeVar2, maxIter,
			// nuround mu
		param[0], nuround, mu, muSq, mup1,
			// cos theta, sin theta
		cos(param[5]), sin(param[5]),
			// parameters from matern.c in geostatsp
			// anisoRatioSq
		(param[4])*(param[4]),
			// varscale
		log(param[2]) - Rf_lgammafn(param[0]) - (param[0]-1)*M_LN2,
			// logxscale
		1.5 * M_LN2 + 0.5 * log(param[0]) - log(param[1]),
			// parameters from bessel temme in gsl
		sinrat, g_1pnu, g_1mnu, g1, g2,
		GSL_DBL_EPSILON /1000, 
		*vclCoords,
		*vclVar));

		// diagonal matrix with nugget + sigmasq

	const double varDiag = param[3] + param[2];
	viennacl::linalg::opencl::matrix_diagonal_assign(*vclVar, varDiag);	

	if( type >= 2 ) {
		// cholesky
		 viennacl::linalg::lu_factorize(*vclVar);
//		vclDofLDL = diag(vclA);
//		for(D=0, D < sizeA1, ++D) {
//			diagOfA[D] = 1.0;
//		}
	}
		// or matrix_diagonal_assign or vector_assign?
		// 3, 4 probably won't be implemente
//	if(type == 5) { // solve for Lt b = X
// not implemented yet
//	}



}

