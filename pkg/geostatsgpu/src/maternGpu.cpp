#include "geostatsgpu.hpp"

#include "materndotcl.hpp"

using namespace Rcpp;
using namespace viennacl;
using namespace viennacl::linalg;


double maternGpuVcl(
	viennacl::matrix<double> &vclVar,
	viennacl::matrix<double> &vclCoords,
	viennacl::vector_base<double> &DofLDL,
	double *param, // range shape variance nugget ratio angleRadians
	const int type,
	viennacl::ocl::kernel &maternKernel
	){


	double logdet=0.0; // the result

	const unsigned int 
	iSizeCoords2=vclCoords.internal_size2(),
	iSizeVar1=vclVar.internal_size1(),
	iSizeVar2=vclVar.internal_size2(),
	sizeVar1=vclVar.size1(),
	Ncell = sizeVar1 * (sizeVar1 - 1)/2,
	maxIter = 1500;


	int nuround = (int) (param[1]+0.5);
	double mu = param[1] - nuround;
	double g_1pnu, g_1mnu, g1, g2;
	const double muSq = mu*mu, 
	varDiag = param[3] + param[2],
	mup1 = mu + 1.0,
	pi_nu = M_PI * mu;
	const double sinrat = (fabs(pi_nu) < GSL_DBL_EPSILON ? 1.0 : pi_nu/sin(pi_nu));

	Rtemme_gamma(&mu, &g_1pnu, &g_1mnu, &g1, &g2);

	// execute kernel

	viennacl::ocl::enqueue(maternKernel(
		Ncell, iSizeCoords2, iSizeVar1, iSizeVar2, maxIter,
			// nuround mu
		param[1], nuround, mu, muSq, mup1,
			// cos theta, sin theta
		cos(param[5]), sin(param[5]),
			// parameters from matern.c in geostatsp
			// anisoRatioSq
		(param[4])*(param[4]),
			// varscale
		log(param[2]) - Rf_lgammafn(param[1]) - (param[1]-1)*M_LN2,
			// logxscale
		1.5 * M_LN2 + 0.5 * log(param[1]) - log(param[0]),
			// parameters from bessel temme in gsl
		sinrat, g_1pnu, g_1mnu, g1, g2,
		GSL_DBL_EPSILON /1000, 
		vclCoords,
		vclVar));

	viennacl::linalg::opencl::matrix_diagonal_assign(vclVar, varDiag);	


	if( type >= 2 ) {
		// cholesky
		logdet = luT<double>(vclVar, DofLDL);
		/*
		viennacl::linalg::lu_factorize(vclVar);
		// try cusolverDnDpotrf instead?

		// pointer to the actual diagonal
		viennacl::vector_base<double> diagOfVar(
			vclVar.handle(), vclVar.size1(), 0, vclVar.internal_size2() + 1);

		// compute log determinant
		DofLDL = element_log(diagOfVar);
		logdet = viennacl::linalg::sum(DofLDL);
// OPERATION_UNARY_LOG_TYPE 	
		//http://viennacl.sourceforge.net/doc/scheduler_8cpp-example.html#a11

		// put the diagonals in D, and 1's on the diagonal of L
		DofLDL = diagOfVar;
		diagOfVar = 1.0;
		*/
	}

	return(logdet);
}


double maternGpuVcl(
	viennacl::matrix<double> &vclVar,
	viennacl::matrix<double> &vclCoords,
	viennacl::vector_base<double> &DofLDL,
	double *param,
	const int type,
	const int ctx_id,
	int max_local_size
	){
	// the context
	viennacl::ocl::context ctx(viennacl::ocl::get_context(ctx_id));

	cl_device_type type_check = ctx.current_device().type();

	// given context but no kernel
	// add kernel to program
	viennacl::ocl::program & my_prog = ctx.add_program(
		maternCLstring,
		"my_kernel");
	// get compiled kernel function
	viennacl::ocl::kernel & maternKernel = my_prog.get_kernel("maternCL");


	// set global work sizes
	const unsigned int 
	sizeVar1=vclVar.size1(),
	sizeVar2=vclVar.size2(),
	Ncell = sizeVar1 * (sizeVar1 - 1)/2;

	if(max_local_size > Ncell) max_local_size = Ncell;

	const double workRatio = Ncell/max_local_size;
	const int workRatioInt = ceil(workRatio);
	int globalSize = workRatioInt*max_local_size;

	// set work sizes
	maternKernel.global_work_size(0, globalSize);
	maternKernel.local_work_size(0, max_local_size);

	double logdet = maternGpuVcl(
		vclVar, vclCoords, DofLDL,
		param, type, maternKernel);
	
	return(logdet);

}

//	SEXP XYR, // solve for Lt b = XY
//	SEXP crossprodR, //bt b


//[[Rcpp::export]]
SEXP cpp_maternGpu(
	Rcpp::S4 varR,
	Rcpp::S4 DofLDLR,
	Rcpp::S4 coordsR,
	Rcpp::NumericVector param, //'range','shape','variance','nugget','anisoRatio','anisoAngleRadians'
	const int type, // 2 cholesky 3 inversecholesky, 4 inverse, 5 solve for b
	int max_local_size) {

	double logdet = 0.0;

	// data
	const bool BisVCL=1;
	const int ctx_id = INTEGER(varR.slot(".context_index"))[0]-1;

	std::shared_ptr<viennacl::matrix<double> > vclVar = getVCLptr<double>(
		varR.slot("address"), BisVCL, ctx_id);


	std::shared_ptr<viennacl::matrix<double> > vclCoords = getVCLptr<double>(
		coordsR.slot("address"), BisVCL, ctx_id);
	// vector to contain the D
	std::shared_ptr<viennacl::vector_base<double> > DofLDL = getVCLVecptr<double>(
		DofLDLR.slot("address"), BisVCL, ctx_id);

	double *param2 = &param[0];

	logdet = maternGpuVcl(
		*vclVar, *vclCoords, *DofLDL, 
		param2,
		type, ctx_id, 
		max_local_size);

	return(Rcpp::wrap(logdet));	
}

