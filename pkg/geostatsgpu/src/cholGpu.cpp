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
	const int colGroupwise,
	int Ncrossprod,
	const int verbose,
	// kernels
	viennacl::ocl::kernel &cholDiag,
	viennacl::ocl::kernel &cholOffDiag,
	viennacl::ocl::kernel &cholOffDiagGroupwise,
	viennacl::ocl::kernel &cholCrossprod,
	viennacl::ocl::kernel &cholSumCrossprod,
	viennacl::ocl::kernel &cholFromCrossprod,
	viennacl::ocl::kernel &sumLog
	){

	double tempDouble; 
	int Dcol, Ncycles, Ncyclesm1, NbeforeLastCycle;
	int Drow, DrowFromZero, Nrows, DcolP1;
	int DcolNpad;
	int rowToStartGroupwise, Dcycle;

	const int NlocalStorage = diagLocal.size() / sizeof(cl_double);

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
	const int NcrossprodI = floor(sqrt(NcrossprodSquared));
	Ncrossprod = std::min(Ncrossprod, NcrossprodI);
	const int NbeginLast = N - Ncrossprod;
	const int NlocalSum = 4;
	const int NrowsInCrossprod = N - NbeginLast;
	const int NrowsInCrossprodSq = 
		NrowsInCrossprod * NrowsInCrossprod;
	const int NbeginLastNpad = NbeginLast * Npad;
	const int Nitemwise = std::min(NbeginLast, colGroupwise);


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

if(verbose){
	Rcpp::Rcout << " N " << N << 
//			" Npad " << Npad<< 
//			" NlocalStorage " << NlocalStorage<< 
//			" NlocalStorageD " << NlocalStorageD<< 
//			" NlocalStorageSqrt " << NlocalStorageSqrt<< 
//			" NcrossprodSquared " << NcrossprodSquared << 
			" Nitemwise " << Nitemwise << 
  	  " NbeginLast " << NbeginLast << 
	    " Ncrossprod " << Ncrossprod << 
			"\n";
}

	for(Dcol=2;Dcol<Nitemwise;Dcol++) {
		DcolNpad = Dcol * Npad;		

		Ncycles = ceil(Dcol / NlocalStorageD);
		Ncyclesm1 = Ncycles-1;

		Nrows = N - Dcol;

		if(verbose) {
			if( (Dcol < 4) | (Dcol < 1005 & Dcol > 1000)) {

			Rcpp::Rcout << "Dcol " << Dcol << 
				" Nrows " << Nrows << 
				" Nglobal " << Nglobal << 
				" floor " << floor(Nrows / Nglobal)<<
				" Ncycles " << Ncycles<<
				" NbeforeLastCycle " << NbeforeLastCycle <<
				"\n";
			}
		}

		
		// diagonals and diagTimesRowOfA
		viennacl::ocl::enqueue(cholDiag(
			A, D,
			diagWorking, 
			diagTimesRowOfA,
			diagLocal, 
			Dcol, Npad, NlocalSum));


		// sum diagWorking to get diag[Dcol]
		tempDouble = A(Dcol,Dcol) - 
			viennacl::linalg::sum(diagWorking);
		D(Dcol) = tempDouble;	

		if(verbose) {
			if( (Dcol < 4) | (Dcol < 1005 & Dcol > 1000)) {
  			Rcpp::Rcout << " diag " << tempDouble <<
	  		" the sum " << viennacl::linalg::sum(diagWorking) <<
		  	"\n";
		  }
		}
		// first cycles through columns
		// cache as many diagonal products as possible
		// in local memory, divide by 1.0
		for(Dcycle = 0; Dcycle < Ncyclesm1; Dcycle++) {
			viennacl::ocl::enqueue(cholOffDiag(
				A, diagTimesRowOfA,
				diagLocal, 1.0,
				Dcol, DcolNpad, N, Npad,
				Dcycle * NlocalStorage,
				(Dcycle +1) * NlocalStorage,
				NlocalStorage,
				NlocalSum));
		}

		// last cycle 
		// fewer items available to cache
		// divide result by D(0)
		if(verbose) {
			if( (Dcol < 5) | (Dcol < 1005 & Dcol > 1000)) {
			Rcpp::Rcout << "Dcol " << Dcol << 
				" startCol " << Dcycle * NlocalStorage << 
				" endCol " << Dcol << 
				" numCol " << Dcol- Dcycle * NlocalStorage <<
				"\n";

		}}


		Dcycle = Ncyclesm1;
		{
			viennacl::ocl::enqueue(cholOffDiag(
				A, diagTimesRowOfA,
				diagLocal, tempDouble,
				Dcol, DcolNpad, N, Npad,
				Dcycle * NlocalStorage,
				Dcol,
				Dcol - Dcycle * NlocalStorage,
				NlocalSum));
		}


	} // end loop through columns
	if(verbose) {
	  	Rcpp::Rcout << "start one row per group\n";
	  }
	int Nverbose = Nitemwise + 10;
	// switch to one row per work group
	for(Dcol = Nitemwise; Dcol < NbeginLast; Dcol++) {

	  if(verbose) {
	    if(Dcol < Nverbose) {
	      Rcpp::Rcout << Dcol;
	    }}
	  
		// diagonals and diagTimesRowOfA
		viennacl::ocl::enqueue(cholDiag(
			A, D,
			diagWorking, 
			diagTimesRowOfA,
			diagLocal, 
			Dcol, Npad, NlocalSum));

	  if(verbose) {
	    if(Dcol < Nverbose) {
	      Rcpp::Rcout << ".";
	    }}
		// sum diagWorking to get diag[Dcol]
		tempDouble = A(Dcol,Dcol) - 
			viennacl::linalg::sum(diagWorking);
		D(Dcol) = tempDouble;	
		if(verbose) {
		  if(Dcol < Nverbose) {
		    Rcpp::Rcout << ".";
		  }}
		viennacl::ocl::enqueue(cholOffDiagGroupwise(
			A, diagTimesRowOfA,
			diagLocal, tempDouble,
			Dcol, Dcol*Npad, N, Npad,
			Ngroups,
			NlocalSum));
		if(verbose) {
		  if(Dcol < Nverbose) {
		    Rcpp::Rcout << ".\n";
		  }}
		
	}

	// last few columns
	// compute the cross product
	// then finish the cholesky
	if(verbose) {
		Rcpp::Rcout << "crossproducts\n";
	}

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

	if(verbose) {
		Rcpp::Rcout << "final columns" << "\n";
	}

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
	const int colGroupwise,
	int Ncrossprod,
	const int verbose,
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
		&cholOffDiagGroupwise = my_prog.get_kernel("cholOffDiagGroupwise"),
		&cholCrossprod = my_prog.get_kernel("cholCrossprod"),
		&cholSumCrossprod = my_prog.get_kernel("cholSumCrossprod"),
		&cholFromCrossprod = my_prog.get_kernel("cholFromCrossprod"),
		&sumLog = my_prog.get_kernel("sumLog");

	// set global work sizes

	cholDiag.global_work_size(0, MCglobal);
	cholDiag.local_work_size(0, MClocal);

	cholOffDiag.global_work_size(0, MCglobal);
	cholOffDiag.local_work_size(0, MClocal);

	cholOffDiagGroupwise.global_work_size(0, MCglobal);
	cholOffDiagGroupwise.local_work_size(0, MClocal);

	sumLog.global_work_size(0, MCglobal);
	sumLog.local_work_size(0, MClocal);

	cholCrossprod.global_work_size(0, MCglobal);
	cholCrossprod.local_work_size(0, MClocal);

	cholSumCrossprod.global_work_size(0, MCglobal);
	cholSumCrossprod.local_work_size(0, MClocal);

	cholFromCrossprod.global_work_size(0, MCglobal);
	cholFromCrossprod.local_work_size(0, MClocal);

if(verbose){

	Rcpp::Rcout <<
		" chol Off global size " << cholOffDiag.global_work_size(0) << 
		" chol Off local size " << cholOffDiag.local_work_size(0) << 
		" local vector size " << localVectorSize << 
	"	\n";
}

	viennacl::ocl::local_mem diagLocal(
		localVectorSize*sizeof(cl_double));


	// one entry per group, for summing
//	viennacl::vector_base<double> Dworking(
//		ceil(MCglobal/MClocal)*sizeof(cl_double));

//	viennacl::vector_base<double> diagTimesRowOfA(
//		x.size1()*sizeof(cl_double));

if(verbose){

	Rcpp::Rcout <<
		" diag local size " << diagLocal.size()/sizeof(cl_double) << 
		"\n";
}

	double logdet = cholGpu(
		x, D, 
		diagWorking, 
		diagTimesRowOfA,
		diagLocal,
		colGroupwise,
		Ncrossprod,
		verbose,
		cholDiag,
		cholOffDiag,
		cholOffDiagGroupwise,
		cholCrossprod,
		cholSumCrossprod,
		cholFromCrossprod,
		sumLog);


	return(logdet);

}


//[[Rcpp::export]]
SEXP cpp_cholGpu(
	Rcpp::S4 xR,
	Rcpp::S4 DR,
	Rcpp::S4 diagWorkingR,
	Rcpp::S4 diagTimesRowOfAR,
	int MCglobal,
	int MClocal,
  	int localStorage,
  	int colGroupwise,
  	int Ncrossprod,
	int verbose,
	std::string kernelR) {

	const int ctx_id = INTEGER(xR.slot(".context_index"))[0]-1;

	const bool BisVCL=1;
	double logdet = -99.9;

	std::shared_ptr<viennacl::matrix<double> > 
		x = getVCLptr<double>(xR.slot("address"), BisVCL, ctx_id);

	std::shared_ptr<viennacl::vector_base<double> > 
		D = getVCLVecptr<double>(
			DR.slot("address"), 
			BisVCL, ctx_id);
	std::shared_ptr<viennacl::vector_base<double> > 
		diagWorking = getVCLVecptr<double>(
			diagWorkingR.slot("address"), 
			BisVCL, ctx_id);
	std::shared_ptr<viennacl::vector_base<double> > 
		diagTimesRowOfA = getVCLVecptr<double>(
			diagTimesRowOfAR.slot("address"), 
			BisVCL, ctx_id);


	logdet = cholGpu(
		*x, 
		*D, 
		*diagWorking,
		*diagTimesRowOfA, 
		kernelR,
		ctx_id, 
		MCglobal,
		MClocal,
		colGroupwise,
		Ncrossprod,
		verbose,
		localStorage);

	return(Rcpp::wrap(logdet));
}