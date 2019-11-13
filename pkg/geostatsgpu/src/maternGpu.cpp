#include "geostatsgpu.hpp"

#include "materndotcl.hpp"
//#include "dmaterndotcl.hpp"

#define GSL_FLT_EPSILON 1.1920928955078125e-07
#define GSL_DBL_EPSILON 2.2204460492503131e-16

using namespace Rcpp;
using namespace viennacl;
using namespace viennacl::linalg;


template <typename T> 
std::string maternClString() {
  return("undefined");}
template <> std::string maternClString<double>(){
  return(maternCLstringDouble);
}
template <> std::string maternClString<float>(){
  return(maternCLstringFloat);
}



template<typename T>
double maternGpuVcl(
	viennacl::matrix<T> &vclVar,
	viennacl::matrix<T> &vclCoords,
	viennacl::vector_base<T> &DofLDL,
	std::vector<double> param, // range shape variance nugget ratio angleRadians
	const int form,
	T epsilon,
	viennacl::ocl::kernel &maternKernel){

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
	const double 
	muSq = mu*mu, 
	varDiag = param[3] + param[2],
	mup1 = mu + 1.0,
	pi_nu = M_PI * mu,
	sinrat = (fabs(pi_nu) < GSL_DBL_EPSILON ? 1.0 : pi_nu/sin(pi_nu));

	Rtemme_gamma(&mu, &g_1pnu, &g_1mnu, &g1, &g2);

//	Rcout << " hi " << vclCoords(0,0) << "  " << vclCoords(0,1) <<
//	  vclCoords(1,0) << "  " << vclCoords(1,1) << " ih \n";

	// execute kernel
	viennacl::ocl::enqueue(maternKernel(
	    Ncell, 
		iSizeCoords2, iSizeVar1, 
		iSizeVar2, maxIter,
		(T) (param[1]), //nu, 
      	(T) (param[5]), // theta
			// parameters from matern.c in geostatsp
			// anisoRatioSq
		(T) ( (param[4])*(param[4]) ),
			// varscale
		(T) (log(param[2]) - Rf_lgammafn(param[1]) - (param[1]-1)*M_LN2),
			// logxscale, includes range parameter
		(T) (1.5 * M_LN2 + 0.5 * log(param[1]) - log(param[0])),
			// parameters from bessel temme in gsl
		(T) (sinrat), 
		(T) (g_1pnu), 
		(T) (g_1mnu), 
		(T) (g1), (T) (g2),
		epsilon,
		vclCoords,
		vclVar));



	viennacl::linalg::opencl::matrix_diagonal_assign(vclVar, (T) (varDiag) );	

	if( form >= 2 ) {
		// cholesky
		logdet = luT<T>(vclVar, DofLDL);
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

// note all objects on the cpu are double
// only the gpu objects are typed
template<typename T> 
double maternGpuVcl(
	viennacl::matrix<T> &vclVar,
	viennacl::matrix<T> &vclCoords,
	viennacl::vector_base<T> &DofLDL,
	std::vector<double> param,
	const int form,
	const int ctx_id,
	std::vector<int> numWorkItems,
	std::vector<int> numLocalItems)
{
	// the context
	viennacl::ocl::context ctx(viennacl::ocl::get_context(ctx_id));

	cl_device_type type_check = ctx.current_device().type();

	std::string maternClStringWithOptions = maternClString<T>();

	viennacl::ocl::program & my_prog = ctx.add_program(maternClStringWithOptions, "my_kernel");
	// get compiled kernel function
	viennacl::ocl::kernel & maternKernel = my_prog.get_kernel("maternCL");

	// set global work sizes
	const unsigned int 
	sizeVar1=vclVar.size1(),
	sizeVar2=vclVar.size2(),
	Ncell = sizeVar1 * (sizeVar1 - 1)/2;

//	if(max_local_size > Ncell) max_local_size = Ncell;

//	const double workRatio = Ncell/max_local_size;
//	const int workRatioInt = ceil(workRatio);
//	int globalSize = workRatioInt*max_local_size;

	// set work sizes
	
//		Rcout << " hi " << numWorkItems[0] << "  " << numWorkItems[1] <<
//		  " " << numWorkItems[2] <<  " ih \n";

	double logdet=0.0;
		
	maternKernel.global_work_size(0, (cl_uint) (numWorkItems[0] ) );//numWorkItems[0]);
	maternKernel.global_work_size(1, (cl_uint) (numWorkItems[1] ) );//numWorkItems[0]);

	maternKernel.local_work_size(0,(cl_uint) (numLocalItems[0]));
	maternKernel.local_work_size(1,(cl_uint) (numLocalItems[1]));
	
	
//	Rcout << " 2hi " << vclCoords(0,0) << "  " << vclCoords(0,1) <<
//	  vclCoords(1,0) << "  " << vclCoords(1,1) << " ih2 \n";

//Rcout << " hi " << maternKernel.global_work_size(0) << "  " << 
//  maternKernel.global_work_size(1) <<
//  " " << maternKernel.global_work_size(2) <<  " ih \n";

	logdet = maternGpuVcl<T>(
		vclVar, vclCoords, DofLDL, 
		param, form, 
		maternClEpsilon<T>(), 
		maternKernel);

	return(logdet);
}

//	SEXP XYR, // solve for Lt b = XY
//	SEXP crossprodR, //bt b
template<typename T> 
SEXP maternGpuVcl(
	Rcpp::S4 varR,
	Rcpp::S4 coordsR,
	Rcpp::S4 DofLDLR,
	Rcpp::NumericVector param, //'range','shape','variance','nugget','anisoRatio','anisoAngleRadians'
	const int form, // 2 cholesky 3 inversecholesky, 4 inverse, 5 solve for b
	Rcpp::IntegerVector numWorkItems,
	Rcpp::IntegerVector numLocalItems) {

	double logdet = 0.0;
	std::vector<double> param2 = Rcpp::as<std::vector<double> >(param);
	std::vector<int> numWorkItemsStd = Rcpp::as<std::vector<int> >(numWorkItems);
	std::vector<int> numLocalItemsStd = Rcpp::as<std::vector<int> >(numLocalItems);
	
	// data
	const bool BisVCL=1;
	const int ctx_id = INTEGER(varR.slot(".context_index"))[0]-1;
	std::shared_ptr<viennacl::matrix<T> > vclVar = getVCLptr<T>(varR.slot("address"), BisVCL, ctx_id);

	// vector to contain the D
	std::shared_ptr<viennacl::vector_base<T> > DofLDL = getVCLVecptr<T>(DofLDLR.slot("address"), BisVCL, ctx_id);

	std::shared_ptr<viennacl::matrix<T> > vclCoords = getVCLptr<T>(coordsR.slot("address"), BisVCL, ctx_id);
	
	logdet = maternGpuVcl<T>(
		*vclVar, *vclCoords, *DofLDL, param2, form, ctx_id,  numWorkItemsStd, numLocalItemsStd);

	return(Rcpp::wrap(logdet));	

	}

//[[Rcpp::export]]
SEXP cpp_maternGpuD(
	Rcpp::S4 varR,
	Rcpp::S4 coordsR,
	Rcpp::S4 DofLDLR,
	Rcpp::NumericVector param, //'range','shape','variance','nugget','anisoRatio','anisoAngleRadians'
	const int form, // 2 cholesky 3 inversecholesky, 4 inverse, 5 solve for b
	Rcpp::IntegerVector numWorkItems,
	Rcpp::IntegerVector numLocalItems) {

	return(maternGpuVcl<double>(
		varR, coordsR, DofLDLR, param, form, numWorkItems, numLocalItems));
}


//[[Rcpp::export]]
SEXP cpp_maternGpuF(
	Rcpp::S4 varR,
	Rcpp::S4 coordsR,
	Rcpp::S4 DofLDLR,
	Rcpp::NumericVector param, //'range','shape','variance','nugget','anisoRatio','anisoAngleRadians'
	const int form, // 2 cholesky 3 inversecholesky, 4 inverse, 5 solve for b
	Rcpp::IntegerVector numWorkItems,
	Rcpp::IntegerVector numLocalItems) {

	return(maternGpuVcl<float>(
		varR, coordsR, DofLDLR, param, form, numWorkItems, numLocalItems));
}




