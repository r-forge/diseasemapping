
#include "geostatsgpu.hpp"

using namespace Rcpp;
using namespace viennacl;


double cholGpu(
	// matrix to cholesky
	viennacl::matrix<double> &A,
	// on exit the diagonal
	viennacl::vector_base<double> &D,
	// temporary storage, of length local_size
	viennacl::vector_base<double> &diagWorking,
	// store diag times rows of A, length N
	viennacl::vector_base<double> &diagTimesRowOfA,
	// store sections of diagTimesRowOfA locally
	viennacl::ocl::local_mem &diagLocal,
	const int NlocalStorage,
	// kernels
	viennacl::ocl::kernel &cholDiag,
	viennacl::ocl::kernel &cholOffDiag,
	viennacl::ocl::kernel &sumLog
	){

	double tempDouble; // the result
	double NlocalStorateD = NlocalStorage;
	int Dcol, Dcolm1, Ncycles, Ncyclesm1, NbeforeLastCycle;
	int err;

	const int 
		Npad=A.internal_size2(),
		N=A.size2();

#ifdef UNDEF
	Rcout << "diagLocalSize " << NlocalStorage <<
		"\n";
#endif

	// first column
	viennacl::vector_base<double> firstColX(
			A.handle(), N, 0, 1);
	tempDouble = firstColX(0);
	D(0) = tempDouble;
	firstColX *= (1 / tempDouble);


// remaining columns
	for(Dcol=1;Dcol<N;Dcol++) {
		
		Dcolm1 = Dcol - 1;
		Ncycles = ceil(Dcol / NlocalStorateD);


		Ncyclesm1 = Ncycles-1;
		NbeforeLastCycle = Ncyclesm1 * NlocalStorage;

#ifdef UNDEF
	Rcout << "Dcol " << Dcol << 
		" NlocalStorage " << NlocalStorateD << 
		" Ncycles " << Ncycles << 
		" Nbefore " <<  NbeforeLastCycle << "\n";
#endif


		
		// diagonals and diagTimesRowOfA
		viennacl::ocl::enqueue(cholDiag(
			A, D,
			diagWorking, 
			diagTimesRowOfA,
			diagLocal, 
			Dcol, Dcolm1, Npad));


		// sum diagWorking to get diag[Dcol]
		tempDouble = A(Dcol,Dcol) - 
			viennacl::linalg::sum(diagWorking);
		D(Dcol) = tempDouble;	

#ifdef UNDEF
		Rcout << "D " << diagWorking(0) << " " 
		<< diagWorking(1) << " " <<
		diagWorking(3) << " " << 
		diagWorking(2) << " " <<
		" A " << A(1,1) << " " << 
		" tempd " << tempDouble << "\n";

		Rcout << "DA " << diagTimesRowOfA(0) << " " 
		<< diagTimesRowOfA(1) << " " <<
		diagTimesRowOfA(3) << " " << 
		diagTimesRowOfA(2) << " " <<
		"\n";
#endif

		viennacl::ocl::enqueue(cholOffDiag(
			A, D,
			diagTimesRowOfA,
			diagLocal, 
			tempDouble,
			Dcol, Dcolm1, Dcol*Npad,
			N, Npad, NlocalStorage,
			Ncycles, Ncyclesm1, 
			NbeforeLastCycle,
			Dcol - NbeforeLastCycle
			));
	}

#ifdef UNDEF
	Rcout << "NlocalSize " << D(10) << 
		" NlocalStorage " << D(11) << 
		" " << D(12) << " " << D(13)<< 
		" Nin " << D(14) << " dd " << D(15) <<
		" Dcy " << D(16) << 
		"\n";
#endif

	viennacl::ocl::enqueue(sumLog(
		D, diagWorking, diagLocal, N
		));
	tempDouble = viennacl::linalg::sum(diagWorking);


	return(tempDouble);

}

double cholGpu(
	viennacl::matrix<double> &x,
	viennacl::vector_base<double> &D,
	std::string kernel,
	const int ctx_id,
	const int MCglobal,
	const int MClocal,
	const int localVectorSize
){

	if(MClocal > localVectorSize) {
		Rcpp::stop("MClocal must be less than localSize");
	}
	// the context
	viennacl::ocl::context ctx(viennacl::ocl::get_context(ctx_id));
	cl_device_type type_check = ctx.current_device().type();

	// given context but no kernel
	// add kernel to program
	viennacl::ocl::program & my_prog = ctx.add_program(
		kernel,
		"my_kernel");

	// get compiled kernel functions
	viennacl::ocl::kernel 
		&cholDiag = my_prog.get_kernel("cholDiag"),
		&cholOffDiag = my_prog.get_kernel("cholOffDiag"),
		&sumLog = my_prog.get_kernel("sumLog");

	// set global work sizes

	cholDiag.global_work_size(0, MCglobal);
	cholDiag.local_work_size(0, MClocal);

	cholOffDiag.global_work_size(0, MCglobal);
	cholOffDiag.local_work_size(0, MClocal);

	sumLog.global_work_size(0, MCglobal);
	sumLog.local_work_size(0, MClocal);

	viennacl::ocl::local_mem Dlocal(
		localVectorSize*sizeof(cl_double));
	viennacl::vector_base<double> Dworking(
		ceil(MCglobal/MClocal)*sizeof(cl_double));
	viennacl::vector_base<double> diagTimesRowOfA(
		x.size1()*sizeof(cl_double));

# ifdef UNDEF
	Rcout <<
		"e Dlocal size " << Dlocal.size() << 
		" global size " << cholOffDiag.global_work_size(0) << 
		" local size " << cholOffDiag.local_work_size(0) << 
		" local vector size" << localVectorSize << "\n";
# endif

	double logdet = cholGpu(
		x, D, 
		Dworking, 
		diagTimesRowOfA,
		Dlocal,
		localVectorSize,
		cholDiag,
		cholOffDiag,
		sumLog);

	return(logdet);

}

//[[Rcpp::export]]
SEXP cpp_cholGpu(
	SEXP xR,
	SEXP DR,
	IntegerVector MCglobal,
	IntegerVector MClocal,
  	IntegerVector localStorage,
	IntegerVector ctx_id,
	CharacterVector kernelR) {

	// data
	const bool BisVCL=1;
	std::shared_ptr<viennacl::matrix<double> > 
		x = getVCLptr<double>(xR, BisVCL, ctx_id[0]);
	std::shared_ptr<viennacl::vector_base<double> > 
		D = getVCLVecptr<double>(DR, BisVCL, ctx_id[0]);

	double logdet = cholGpu(
		*x, 
		*D, 
		Rcpp::as< std::string >(kernelR(0)),
		ctx_id[0], 
		MCglobal[0],
		MClocal[0],
		localStorage[0]);

	return(Rcpp::wrap(logdet));
}