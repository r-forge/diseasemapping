
// https://en.wikipedia.org/wiki/Cholesky_decomposition#LDL_decomposition_2

// 1d with 64 work items, 16 groups of 4
// 2d with 64 items, 16 by 1 workgroups, each is 2 by 2

__kernel void cholDiag(
	__global double *A, 
	__global double *diag,
	__global double *diagWorking,
	__global double *diagTimesRowOfA,
	__local double *diagLocal,
	const int Dcol,
	const int Npad) {

	const int Nsize = get_global_size(0);
	const int Dlocal1 = get_local_id(0);
	const int Dlocal2 = get_local_id(1);
	const int NlocalSize1 = get_local_size(0);
	const int NlocalSize2 = get_local_size(1);

	const int DlocalUnique = Dlocal1 + 
		NlocalSize1 * Dlocal2;
	const int Dstart = Dlocal2 * NlocalSize1;
	const int Dend = Dstart + NlocalSize1;

	int Dk;
	double DL, Ddiag, diagDkTimesDl;

	Ddiag = 0.0;
	for(Dk = get_global_id(0); Dk < Dcol; Dk += Nsize) {
		DL = A[Dcol + Dk * Npad];
		diagDkTimesDl = diag[Dk] * DL;

		diagTimesRowOfA[Dk] = diagDkTimesDl;

		Ddiag += diagDkTimesDl * DL;
	}
	diagLocal[DlocalUnique] = Ddiag;
	// must be added to get D[Dcol]
	barrier(CLK_LOCAL_MEM_FENCE);

	// add up diagLocal over the 'columns'
	// of the work groups
	if(Dlocal1==0) {
		Ddiag = 0.0;
		for(Dk = Dstart; Dk < Dend; Dk++) {
			Ddiag += diagLocal[Dk];
		}
		diagLocal[Dlocal1] = Ddiag;
	}
	barrier(CLK_LOCAL_MEM_FENCE);

	// now add the column sums
	if(DlocalUnique==0) {
		Ddiag = 0.0;
		for(Dk = 0; Dk < NlocalSize2; Dk++) {
			Ddiag += diagLocal[Dk]
		}
	}
	diagWorking[get_group_id(0)] = Ddiag;
	// must be added across work groups(0)
	// to get D[Dcol]
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
	const int Dlocal1 = get_local_id(0);
	const int Dlocal2 = get_local_id(1);
	const int NlocalSize1 = get_local_size(0);
	const int NlocalSize2 = get_local_size(1);
	const int DlocalUnique = Dlocal1 + 
		NlocalSize1 * Dlocal2;

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
