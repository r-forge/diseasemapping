
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
	const int Dcolm1,
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
	const int Dcolm1,
	const int DcolNpad,
	const int N, 
	const int Npad,
	const int NlocalStorage,
	const int Ncycles, // ceil(Dcolm1 / NlocalStorage)
	const int Ncyclesm1,
	const int NbeforeLastCycle, // = Ncyclesm1 * NlocalStorage;
	// number of elements in last cycle
	// Dcol - NbeforeLastCycle;
	const int NinLastCycle) {

	const int Dlocal = get_local_id(0);

	const int Dglobal = get_global_id(0);
	const int Nsize = get_global_size(0);
	const int NlocalSize = get_local_size(0);

	int Dcycle, DcycleNlocalStorage;
	int Drow, Dk, DrowDk;
	double DL;
	// Dcycle through groups of NlocalStorage columns
	for(Dcycle=0;Dcycle<Ncyclesm1;Dcycle++) {
//	for(Dcycle=0;Dcycle<1;Dcycle++) {
		DcycleNlocalStorage = Dcycle * NlocalStorage;
		// copy over some of D L rows

		Dk = Dlocal;
		diagLocal[Dk] = diagTimesRowOfA[
				Dk + DcycleNlocalStorage];

		for(Dk = Dlocal; 
			Dk < NlocalStorage; 
			Dk += NlocalSize) {
	
			diagLocal[Dk] = diagTimesRowOfA[
				Dk + DcycleNlocalStorage];

		}

		barrier(CLK_LOCAL_MEM_FENCE);

		// loop through rows
		for(Drow = Dcol; Drow < N; Drow ++) {

			DL = A[Drow + DcolNpad];
			for(Dk = 0; Dk < NlocalStorage; Dk++) {
				DL -= A[
					Drow + (Dk + DcycleNlocalStorage) * Npad
					] * diagLocal[Dk];
			}
			A[Drow + DcolNpad] = DL;
		}

	}

	// the final group of columns
	// has fewer than NlocalStorage elements
	DcycleNlocalStorage = Ncyclesm1 * NlocalStorage;

	for(Dk = Dlocal; 
		Dk < NinLastCycle; 
		Dk += NlocalSize) {
	
		diagLocal[Dk] = diagTimesRowOfA[
				Dk + DcycleNlocalStorage];
	
	}
	barrier(CLK_LOCAL_MEM_FENCE);


	// loop through rows
	for(Drow = Dcol; Drow < N; Drow ++) {

		DL = A[Drow + DcolNpad];
		for(Dk = 0; Dk < NinLastCycle; Dk++) {
			DL -= A[
				Drow + (Dk + DcycleNlocalStorage) * Npad
				] * diagLocal[Dk];
		}
		A[Drow + DcolNpad] = DL / diagDcol;
	}
#ifdef UNDEF
		if(Dglobal==0) {
			diag[10] = NlocalSize;
			diag[11] = NlocalStorage;
			diag[12] = Ncycles;
			diag[13] = Ncyclesm1;
			diag[14] = NinLastCycle;
			diag[15] = diagDcol;
			diag[16] = DcycleNlocalStorage;
		}
#endif
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
