/*
	TO DO:
	- update the R interface
	- make the update one row per work group for the last few rows
	- fully implement the 'type' option in matern
	- clean up unused variables
*/


#include "geostatsgpu.hpp"

using namespace viennacl;


double cholGpu(
	// matrix to cholesky
	viennacl::matrix<double> &A,
	// on exit the diagonal
	viennacl::vector_base<double> &D,
	// temporary storage, of length number of groups
	viennacl::vector_base<double> &diagWorking,
	// store diag times rows of A, length N
	viennacl::vector_base<double> &diagTimesRowOfA,
	// store sections of diagTimesRowOfA locally
	viennacl::ocl::local_mem &diagLocal,
	const int NlocalStorage,
	// kernels
	viennacl::ocl::kernel &cholDiag,
	viennacl::ocl::kernel &cholOffDiag,
	viennacl::ocl::kernel &cholCrossprod,
	viennacl::ocl::kernel &cholSumCrossprod,
	viennacl::ocl::kernel &cholFromCrossprod,
	viennacl::ocl::kernel &sumLog
	){

	double tempDouble; 
	int Dcol, Ncycles, Ncyclesm1, NbeforeLastCycle;
	int Drow, DrowFromZero, Nrows, DcolP1;
	int DcolNpad;

	const int 
		Npad=A.internal_size2(),
		N=A.size2();
	const int Nm1 = N - 1, Nm2 = N-2;

	// how many cross products can be stored in local memory
	const int Nlocal = cholOffDiag.local_work_size();
	const int Nglobal = cholOffDiag.global_work_size(0);			

	const int Ngroups = //cholCrossprod.get_num_groups(0);
		cholCrossprod.global_work_size(0) / 
				cholCrossprod.local_work_size(0);

	const double NlocalStorageD = NlocalStorage;
	const double NlocalStorageSqrt = sqrt(NlocalStorage);

	const double NcrossprodSquared = NlocalStorageD/Nlocal;
	const int Ncrossprod = floor(sqrt(NcrossprodSquared));
	const int NbeginLast = N - Ncrossprod;
	const int NlocalSum = 4;
	const int NrowsInCrossprod = N - NbeginLast;
	const int NrowsInCrossprodSq = 
		NrowsInCrossprod * NrowsInCrossprod;
	const int NbeginLastNpad = NbeginLast * Npad;

	int rowToStartGroupwise = N;

//	viennacl::vector_base<double> oneColX(
//			A.handle(), Nm1, 
//			0, 1);

	viennacl::range rowsFrom1(1, N), rowsFrom2(2,N);
	viennacl::range col1(0, 1), col2(1,2);

	// first column
	viennacl::matrix_range<viennacl::matrix<double>>  
		oneColX(A, col1, rowsFrom1);


	tempDouble = A(0,0);
	D(0) = tempDouble;

	oneColX /= tempDouble;

	// second column
	viennacl::matrix_range<viennacl::matrix<double>>  
		twoColX(A, col2, rowsFrom2);

	// first column, from row 3
	viennacl::matrix_range<viennacl::matrix<double>>
		oneColX2(A, col1, rowsFrom2);

	tempDouble = A(0,1);
	tempDouble = A(1,1) - 
		tempDouble*tempDouble * D(0);
	D(1) = tempDouble;

	twoColX = ( -D(0) * A(0,1) / tempDouble ) * 
//		viennacl::matrix_range<viennacl::matrix<double>>(
//			A, col2, rowsFrom2) + 
	oneColX2 + 
			twoColX / tempDouble; 

	// remaining columns
	// until the number of columns remaining
	// is small enough to store in local memory

#ifdef DEBUG
	Rcpp::Rcout << " N " << N << 
			" Npad " << Npad<< 
			" NlocalStorage " << NlocalStorage<< 
			" NlocalStorageD " << NlocalStorageD<< 
			" NlocalStorageSqrt " << NlocalStorageSqrt<< 
			"\n" <<
			" NbeginLast " << NbeginLast << 
			" NcrossprodSquared " << NcrossprodSquared << 
			" Ncrossprod " << Ncrossprod << 
			"\n";
#endif

	for(Dcol=2;Dcol<NbeginLast;Dcol++) {
		
//		Dcolm1 = Dcol - 1;
		Ncycles = ceil(Dcol / NlocalStorageD);


		Ncyclesm1 = Ncycles-1;
		NbeforeLastCycle = Ncyclesm1 * NlocalStorage;

		Nrows = N - Dcol;
		// row to start one row per workgroup
		rowToStartGroupwise = Dcol + floor(Nrows/Nglobal);
		// unless there are fewer columns to sum over than
		// local work items, or only a small number
		// of work items would be unused if
		// one row per work item were maintained
		if(Dcol < Nlocal | 
			(N - rowToStartGroupwise) < 4
			) {

			rowToStartGroupwise = N;
		}
		
		// diagonals and diagTimesRowOfA
		viennacl::ocl::enqueue(cholDiag(
			A, D,
			diagWorking, 
			diagTimesRowOfA,
			diagLocal, 
			Dcol, Npad));


		// sum diagWorking to get diag[Dcol]
		tempDouble = A(Dcol,Dcol) - 
			viennacl::linalg::sum(diagWorking);
		D(Dcol) = tempDouble;	


		viennacl::ocl::enqueue(cholOffDiag(
			A, D,
			diagTimesRowOfA,
			diagLocal, 
			tempDouble,
			Dcol, Dcol*Npad,
			N, Npad, NlocalStorage,
			Ncyclesm1, 
			NbeforeLastCycle,
			Dcol - NbeforeLastCycle,
			rowToStartGroupwise, 
			NlocalSum));
	}


	// last few columns
	// compute the cross product
	// then finish the cholesky


	// compute the cross products
	viennacl::ocl::enqueue(cholCrossprod(
		A, D, diagTimesRowOfA, diagLocal,
		NbeginLast, NbeginLast * Npad,
		NrowsInCrossprod, NrowsInCrossprodSq,
		Npad, NlocalSum));

	// sum up the cross products from 
	// each work group
	viennacl::ocl::enqueue(cholSumCrossprod(
		diagTimesRowOfA,
 		NrowsInCrossprod, 
		NrowsInCrossprodSq,
		Ngroups));


	for(Dcol = NbeginLast; Dcol < N; Dcol++) {

		DcolP1 = Dcol + 1;
		Nrows = N - Dcol;

//		DcolNpad = Dcol * Npad;



		viennacl::ocl::enqueue(cholFromCrossprod(
			A, D, diagTimesRowOfA, diagLocal,
			Dcol, 
			NbeginLast,
			Dcol - NbeginLast, // colSinceCrossprod
			Dcol + NbeginLastNpad,
			Dcol*Npad + Dcol, // AforUpdateStart
			Ncrossprod,
			(Ncrossprod +1) * (Dcol - NbeginLast),
			Nrows, 
			Npad
			));

		// divide A[,Dcol] by diagonal
		viennacl::matrix_range<viennacl::matrix<double>> theColHere(
			A,
			viennacl::range(Dcol,DcolP1),
			viennacl::range(DcolP1, N));
		theColHere  *= (1/D[Dcol]);
	}

	viennacl::linalg::opencl::matrix_diagonal_assign(A, 1.0);

	viennacl::ocl::enqueue(sumLog(
		D, diagWorking, diagLocal, N,
		NlocalSum
		));

	return(viennacl::linalg::sum(diagWorking));
}


double cholGpu(
	viennacl::matrix<double> &x,
	viennacl::vector_base<double> &D,
	viennacl::vector_base<double> &diagWorking,
	viennacl::vector_base<double> &diagTimesRowOfA,
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
		&cholCrossprod = my_prog.get_kernel("cholCrossprod"),
		&cholSumCrossprod = my_prog.get_kernel("cholSumCrossprod"),
		&cholFromCrossprod = my_prog.get_kernel("cholFromCrossprod"),
		&sumLog = my_prog.get_kernel("sumLog");

	// set global work sizes

	cholDiag.global_work_size(0, MCglobal);
	cholDiag.local_work_size(0, MClocal);

	cholOffDiag.global_work_size(0, MCglobal);
	cholOffDiag.local_work_size(0, MClocal);

	sumLog.global_work_size(0, MCglobal);
	sumLog.local_work_size(0, MClocal);

	cholCrossprod.global_work_size(0, MCglobal);
	cholCrossprod.local_work_size(0, MClocal);

	cholSumCrossprod.global_work_size(0, MCglobal);
	cholSumCrossprod.local_work_size(0, MClocal);

	cholFromCrossprod.global_work_size(0, MCglobal);
	cholFromCrossprod.local_work_size(0, MClocal);


	viennacl::ocl::local_mem diagLocal(
		localVectorSize*sizeof(cl_double));

	// one entry per group, for summing
//	viennacl::vector_base<double> Dworking(
//		ceil(MCglobal/MClocal)*sizeof(cl_double));

//	viennacl::vector_base<double> diagTimesRowOfA(
//		x.size1()*sizeof(cl_double));

# ifdef DEBUG
	Rcpp::Rcout <<
		"e Dlocal size " << diagLocal.size() << 
		" global size " << cholOffDiag.global_work_size(0) << 
		" local size " << cholOffDiag.local_work_size(0) << 
		" local vector size" << localVectorSize << "\n";
# endif

	double logdet = cholGpu(
		x, D, 
		diagWorking, 
		diagTimesRowOfA,
		diagLocal,
		localVectorSize,
		cholDiag,
		cholOffDiag,
		cholCrossprod,
		cholSumCrossprod,
		cholFromCrossprod,
		sumLog);

	return(logdet);

}


//[[Rcpp::export]]
SEXP cpp_cholGpu(
	SEXP xR,
	SEXP DR,
	SEXP diagWorkingR,
	SEXP diagTimesRowOfAR,
	int MCglobal,
	int MClocal,
  	int localStorage,
	int ctx_id,
	std::string kernelR) {


	const bool BisVCL=1;
	std::shared_ptr<viennacl::matrix<double> > 
		x = getVCLptr<double>(xR, BisVCL, ctx_id);
	std::shared_ptr<viennacl::vector_base<double> > 
		D = getVCLVecptr<double>(
			DR, 
			BisVCL, ctx_id);
	std::shared_ptr<viennacl::vector_base<double> > 
		diagWorking = getVCLVecptr<double>(
			diagWorkingR, 
			BisVCL, ctx_id);
	std::shared_ptr<viennacl::vector_base<double> > 
		diagTimesRowOfA = getVCLVecptr<double>(
			diagTimesRowOfAR, 
			BisVCL, ctx_id);

	double logdet = cholGpu(
		*x, 
		*D, 
		*diagWorking,
		*diagTimesRowOfA, 
		kernelR,
		ctx_id, 
		MCglobal,
		MClocal,
		localStorage);

	return(Rcpp::wrap(logdet));
}