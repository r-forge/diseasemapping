

#include <Rmath.h>


//#define VIENNACL_WITH_CUDA
//#define VIENNACL_WITH_OPENCL

#define GSL_DBL_EPSILON 2.2204460492503131e-16
#define GSL_SQRT_DBL_MAX 1.3407807929942596e+154

// make link to /usr/lib/x86_64-linux-gnu/libOpenCL.so

#include <Rcpp.h>
//#include <RcppEigen.h>


//#include "viennacl/backend/cuda.hpp"
#include "gpuR/getVCLptr.hpp"
#include "gpuR/dynVCLMat.hpp"
#include "gpuR/dynVCLVec.hpp"
#include "viennacl/ocl/backend.hpp"
#include "viennacl/linalg/sum.hpp"
#include "viennacl/linalg/lu.hpp"


extern "C" void Rtemme_gamma(double *nu, double * g_1pnu, double * g_1mnu, double *g1, double *g2);

using namespace Rcpp;
using namespace viennacl;
using namespace viennacl::linalg;

//[[Rcpp::export]]
SEXP cpp_maternGpu(
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

	double logdet=0.0;


	std::string my_kernel = as<std::string>(sourceCode_);
	viennacl::ocl::context ctx(viennacl::ocl::get_context(ctx_id));

	// add kernel to program
	viennacl::ocl::program & my_prog = ctx.add_program(my_kernel, "my_kernel");

	// get compiled kernel function
	viennacl::ocl::kernel & maternKernel = my_prog.get_kernel("maternCL");

	cl_device_type type_check = ctx.current_device().type();

	// data
	const bool BisVCL=1;
	unsigned int D;

	std::shared_ptr<viennacl::matrix<double> > vclVar2 = getVCLptr<double>(varR, BisVCL, ctx_id);
	std::shared_ptr<viennacl::matrix<double> > vclCoords = getVCLptr<double>(coordsR, BisVCL, ctx_id);


	const unsigned int 
	sizeCoords1=vclCoords->size1(),
	sizeCoords2=vclCoords->size2(),
	sizeVar1=vclVar2->size1(),
	sizeVar2=vclVar2->size2(),
	iSizeCoords1=vclCoords->internal_size1(),
	iSizeCoords2=vclCoords->internal_size2(),
	iSizeVar1=vclVar2->internal_size1(),
	iSizeVar2=vclVar2->internal_size2(),
	Ncell = sizeVar1 * (sizeVar1 - 1)/2;

	if(max_local_size > Ncell) max_local_size = Ncell;

	// set global work sizes
	const double workRatio = Ncell/max_local_size;
	const int workRatioInt = ceil(workRatio);
	int globalSize = workRatioInt*max_local_size;


	// set work sizes
	maternKernel.global_work_size(0, globalSize);
	maternKernel.local_work_size(0, max_local_size);

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
		*vclVar2));

		// diagonal matrix with nugget + sigmasq

	const double varDiag = param[3] + param[2];
	viennacl::linalg::opencl::matrix_diagonal_assign(*vclVar2, varDiag);	

	if( type >= 2 ) {
		// cholesky
		viennacl::linalg::lu_factorize(*vclVar2);
		// try cusolverDnDpotrf instead?

		// pointer to the actual diagonal
		viennacl::vector_base<double> diagOfVar( (*vclVar2).handle(), (*vclVar2).size1(), 0, (*vclVar2).internal_size2() + 1);
		// vector to contain the D
		std::shared_ptr<viennacl::vector_base<double> > DofLDL = getVCLVecptr<double>(DofLDLR, BisVCL, ctx_id);

		// compute log determinant
		*DofLDL = element_log(diagOfVar);
		logdet = viennacl::linalg::sum(*DofLDL);
// OPERATION_UNARY_LOG_TYPE 	
		//http://viennacl.sourceforge.net/doc/scheduler_8cpp-example.html#a11

		// put the diagonals in D, and 1's on the diagonal of L
		*DofLDL = diagOfVar;
		diagOfVar = 1.0;

	}
	// 3 or 4 or 5 use
	// viennacl::linalg::inplace_solve	
//	if(type == 5) { // solve for Lt b = X
// not implemented yet
//	}
	return(Rcpp::wrap(logdet));	
}

