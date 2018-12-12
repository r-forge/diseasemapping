
// https://en.wikipedia.org/wiki/Cholesky_decomposition#LDL_decomposition_2

// on exit, diagWorking[i] is 
// sum{A[j]^2 * D[j] ; j in Dlocal + Dk * Nsize}

// diagTimesRowOfA =  diag * A[Dcol, ]

__kernel void cholDiag(
	__global double *A, 
	__global double *diag,
	__global double *diagWorking,
	__global double *diagTimesRowOfA,
	__local double *diagLocal,
	const int Dcol,
	const int Npad) {

	const int Dglobal = get_global_id(0);
	const int Dlocal = get_local_id(0);

	const int Nsize = get_global_size(0);
	const int NlocalSize = get_local_size(0);
	int Dk;
	double DL, Ddiag, diagDkTimesDl;

	Ddiag = 0.0;
	for(Dk = Dglobal; Dk < Dcol; Dk += Nsize) {
		DL = A[Dcol + Dk * Npad];
		diagDkTimesDl = diag[Dk] * DL;

		diagTimesRowOfA[Dk] = diagDkTimesDl;

		Ddiag += diagDkTimesDl * DL;
	}
	diagLocal[Dlocal] = Ddiag;
	// must be added to get D[Dcol]
	barrier(CLK_LOCAL_MEM_FENCE);

	if(Dlocal==0) {
		Ddiag = 0.0;
		for(Dk = 0; Dk < NlocalSize; Dk++) {
			Ddiag += diagLocal[Dk];
		}
		diagWorking[get_group_id(0)] = Ddiag;
	// must be added across work groups to get D[Dcol]
	}

}

__kernel void cholOffDiag(
	__global double *A, 
	__global double *diag,
	__global double *diagTimesRowOfA,
	__local double *diagLocal,
	const double diagDcol,
	const int Dcol,
	const int DcolNpad,
	const int N, 
	const int Npad,
	const int NlocalStorage,
	const int Ncyclesm1,
	const int NbeforeLastCycle, // = Ncyclesm1 * NlocalStorage;
	// number of elements in last cycle
	// Dcol - NbeforeLastCycle;
	const int NinLastCycle) {


	const int Nsize = get_global_size(0);
	const int NlocalSize = get_local_size(0);

	const int Dlocal = get_local_id(0);
	const int DrowStart = Dcol+1+get_global_id(0);

	int Dcycle, DcycleNlocalStorage, Drow, Dk;
	double DL;
	// Dcycle through groups of NlocalStorage columns
	for(Dcycle=0;Dcycle<Ncyclesm1;Dcycle++) {

		DcycleNlocalStorage = Dcycle * NlocalStorage;
		// copy over some of D L rows

		for(Dk = Dlocal; Dk < NlocalStorage; Dk += NlocalSize) {
	
			diagLocal[Dk] = diagTimesRowOfA[
				Dk + DcycleNlocalStorage];

		}
		// synchronize all the work items
		// so that all of diagLocal is ready
		barrier(CLK_LOCAL_MEM_FENCE);

		// loop through rows
		// do this cycle's bit of the cholesky
		for(Drow = DrowStart; Drow < N; Drow+= Nsize) {

			DL = A[Drow + DcolNpad];
			for(Dk = 0; Dk < NlocalStorage; Dk++) {
				DL -= A[
					Drow + (Dk + DcycleNlocalStorage) * Npad
					] * diagLocal[Dk];
			}
			A[Drow + DcolNpad] = DL;
		}
	}
	barrier(CLK_LOCAL_MEM_FENCE);

	// the final group of columns
	// has fewer than NlocalStorage elements
//	DcycleNlocalStorage = = NbeforeLastCycle = Ncyclesm1 * NlocalStorage;

	for(Dk = Dlocal; Dk < NinLastCycle; Dk += NlocalSize) {
	
		diagLocal[Dk] = diagTimesRowOfA[
				Dk + NbeforeLastCycle];//DcycleNlocalStorage];
	
	}
	barrier(CLK_LOCAL_MEM_FENCE);

	for(Drow = DrowStart; Drow < N; Drow+= Nsize) {

		DL = A[Drow + DcolNpad];
		for(Dk = 0; Dk < NinLastCycle; Dk++) {
			DL -= A[
				Drow + (Dk + NbeforeLastCycle) * Npad
				] * diagLocal[Dk];
		}
		A[Drow + DcolNpad] = DL / diagDcol;
	}
}

__kernel void sumLog(
	__global double *D, 
	__global double *diagWorking,
	__local double *diagLocal,
	const int N) {

	const int Dlocal = get_local_id(0);
	const int Nsize = get_global_size(0);
	const int NlocalSize = get_local_size(0);

	int Dk;
	double Ddiag = 0.0;

	for(Dk = get_global_id(0); Dk < N; Dk += Nsize) {
		Ddiag += log(D[Dk]);
	}	
	diagLocal[Dlocal] = Ddiag;
	barrier(CLK_LOCAL_MEM_FENCE);

	// add the local storage
	if(Dlocal==0) {
		Ddiag = 0.0;
		for(Dk = 0; Dk < NlocalSize; Dk++) {
			Ddiag += diagLocal[Dk];
		}
	diagWorking[get_group_id(0)] = Ddiag;
	}

}


// crossprod_ij contains sum_k L_kj L_ij D_k
// crossprodLocal must contain over Nlocal * Nrows^2 
// entries
__kernel void cholCrossprod(
	__global double *A, 
	__global double *diag,
	__global double *crossprod,
	__local double *crossprodLocal,
	const int Dcol,
	const int DcolNpad,
	const int Nrows, // N - Dcol
	const int NrowsSq,
	const int NlocalSum// typicaly 4, 
	) {

	const int Dlocal = get_local_id(0);
	const int Nlocal = get_local_size(0);
	const int Dglobal = get_global_id(0);
	const int Nglobal = get_global_size(0);

// location in crossprodLocal to store the results
	const int crossprodIndex = Dlocal*NrowsSq;

	int Drow1, Drow2, Dk, Drow1Nrows, DkNpad, DlocalNsq;
	double ADhere, Dhere;

	// Dk loops through columns of A
	// do the first column, 
	// to initialize crossprodLocal
	Dk = Dglobal;
	Dhere = diag[Dk];
	DkNpad = Dk * Npad + Dcol;

	// Drow1 loops through rows from Dcol
	for(Drow1 = 0; Drow1 < Nrows; Drow1++) {
		Drow1Nrows = crossprodIndex + Drow1 * Nrows;
		ADhere = A[DkNpad + Drow1] * Dhere;

		for(Drow2 = Drow1; Drow2 < Nrows; Drow2++){
			crossprodLocal[Drow1Nrows + Drow2] =
				ADhere * A[DkNpad + Drow2]; 
		}
	}

	// now loop through the rest of the columns
	for(Dk = Dglobal+ Nglobal; Dk < Dcol; Dk += Nglobal) {
		Dhere = diag[Dk];
		DkNpad = Dk * Npad + Dcol;
		for(Drow1 = 0; Drow1 < Nrows; Drow1++) {
			Drow1Nrows = crossprodIndex + Drow1 * Nrows;
			ADhere = A[DkNpad + Drow1] * Dhere;
			for(Drow2 = Drow1; Drow2 < Nrows; Drow2++){
				crossprodLocal[Drow1Nrows + Drow2] +=
					ADhere * A[DkNpad + Drow2]; 
			}
	}
	barrier(CLK_LOCAL_MEM_FENCE);

	// sum over the local work items
	if(Dlocal < NlocalSum) {
		DlocalNsq = Dlocal * NrowsSq;
		for(Dk = Dlocal + NlocalSum; Dk < Nlocal; 
			Dk += NlocalSum) {
			DkNpad = Dk * NrowsSq;

			for(Drow1 = 0; Drow1 < Nrows; Drow1++) {
			Drow1Nrows = Drow1 * NrowsSq;	
			for(Drow2 = Drow1; Drow2 < Nrows; Drow2++){
				crossprodLocal[
					DlocalNsq + Drow1Nrows+ Drow2
				] += crossprodLocal[
					DlocalNsq +  
					DkNpad + 
					Drow1Nrows +
					Drow2
				];
			}
			}
		}
	}
	barrier(CLK_LOCAL_MEM_FENCE);

	// sum over the local work items
	if(Dlocal == 0 ) {
		for(Dk = 1; Dk < NlocalSum; Dk ++) {
			DkNpad = Dk * NrowsSq;

			for(Drow1 = 0; Drow1 < Nrows; Drow1++) {
			Drow1Nrows = Drow1 * NrowsSq;	
			for(Drow2 = Drow1; Drow2 < Nrows; Drow2++){
				crossprodLocal[
					Drow1Nrows + Drow2
				] += crossprodLocal[
					DkNpad + 
					Drow1Nrows +
					Drow2
				];
			}
			}
		}
		// copy to global memory
		DkNpad = get_group_id(0) * NrowsSq;
	
		for(Drow1 = 0; Drow1 < Nrows; Drow1++) {
			Drow1Nrows = Drow1 * NrowsSq;	
			for(Drow2 = Drow1; Drow2 < Nrows; Drow2++){
				crossprod[
					DkNpad + Drow1Nrows + Drow2
				] += crossprodLocal[
					Drow1Nrows +
					Drow2
				];
			}
		}
	} // end Dlocal == 0

}

__kernel void cholSumCrossprod(
	__global double *crossprod,
	const int Nrows,
	const int NrowsSq,
	const int NtoAdd // number of groups in cholCrossprod
	) {

	const int Dlocal = get_local_id(0);
	const int Dgroup = get_group_id(0);
	const int Ngroups = get_group_size(0);
	const int Nlocal = get_local_size(0);

	int Drow1, Drow2, Dindex, Dsum;
	double theSum;

	for(Drow1 = Dgroup; Drow < Nrows; Drow += Ngroups) {
		for(Drow2 = Drow1 + Dlocal; Drow2 < Nrows; 
			Drow2 += Nlocal) {
			
			Dindex = Drow1 * Nrows + Drow2;
			theSum = 0.0;
			
			for(Dsum = 0; Dsum < NtoAdd; Dsum++) {
				theSum += crossprod[
					Dsum * NrowsSq + Dindex];
			}
			
			crossprod[Dindex] = theSum;
		}
	}
}

__kernel void cholFromCrossprod(
	__global double *A, 
	__global double *diag,
	__global double *crossprod,
	__local double *diagLocal,

	// column of A to update
	const int Dcol,

	// column up to which the cross product was computed
	const int colAtCrossprod,

	// Dcol - colAtCrossprod
	const int colSinceCrossprod,

	// index for A[Dcol, colSinceCrossprod]
	// = Dcol + colSinceCrossprod * Npad;
	const int AforCrossprodStartJ,

	// index for A[Dcol+1, Dcol]
	const int AforUpdateStart,

	const int NrowsCrossprod, // N - colAtCrossprod
	const int Nrows, // N - Dcol
	const int Npad
	) {

	const int Nglobal = get_global_size(0);
	const int Dglobal = get_global_id(0);
	const int Nlocal = get_local_size(0);

	double Dcrossprod;
	int Drow, Dk, AforCrossprodStartI;

	for(Dk = get_local_id(0); Dk < colSinceCrossprod;
		Dk+=Nlocal) {

		Dcrossprod = A[
				AforCrossprodStartJ + Dk * Npad
				] * diag[colAtCrossprod + Dk];
		diagLocal[Dk] = Dcrossprod;
	}
	barrier(CLK_LOCAL_MEM_FENCE);

	// compute the diagonal
	if(Dlocal == 0) {
		diagHere = 0.0;
		for(Dk = 0; Dk < colSinceCrossprod; Dk++) {
			diagHere += diagLocal[Dk];
		}
		diagLocal[colSinceCrossprod + 1] = diagHere;
	}
	barrier(CLK_LOCAL_MEM_FENCE);

	diagHere = diagLocal[colSinceCrossprod + 1];

	if(Dglobal == 0) {
		diag[Dcol] = diagHere;
	}


	for(Drow = Dglobal+1; Drow < Nrows; 
		Drow += Nglobal) {

		// update the crossprod
		AforCrossprodStartI = 
			AforCrossprodStartJ + Drow;
		Dcrossprod = crossprod[
			colSinceCrossprod + Drow];
		for(Dk = 0; Dk < colSinceCrossprod; Dk++) {
			Dcrossprod += A[
				AforCrossprodStartI + Dk * Npad
				] * diagLocal[Dk];
		}

		Dk = AforUpdateStart + Drow;
		A[Dk] = (A[Dk] - Dcrossprod) / diagHere;
	}

}
