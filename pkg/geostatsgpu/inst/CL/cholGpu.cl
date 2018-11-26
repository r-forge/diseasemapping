
// https://en.wikipedia.org/wiki/Cholesky_decomposition#LDL_decomposition_2

__kernel void cholGpu(
	__global double *A, 
	__global double *diag,
	__local double *diagLocal,
	const unsigned int Dcol,
	const unsigned int N, 
	const unsigned int Npad)
{


	const unsigned int Dgroup = get_group_id(0);
	const unsigned int Dlocal = get_local_id(0);
	const unsigned int Dglobal = get_global_id(0);
	const unsigned int Nsize = get_global_size(0);
	const unsigned int NlocalSize = get_local_size(0);
	const unsigned int Ngroups = Nsize / NlocalSize;

	const unsigned int Dcolm1 = Dcol - 1;
	const unsigned int DcolNpad = Dcol * Npad;

	int Drow, Dk, DrowDk;
	const int DcolS = Dcol;

	__local double Ddiag;
	double DL;

	for(Dk = Dlocal; Dk < Dcol; Dk += NlocalSize) {
		DL = A[Dcol + Dk * Npad];
		DL *= DL;
		diagLocal[Dk] = diag[Dk] * DL;
	}

	barrier(CLK_LOCAL_MEM_FENCE);
// add 'em up to get diagonal
	if(Dlocal==0) {
		Ddiag = 0.0;
		for(Dk = 0; Dk < Dcol; Dk++) {
			Ddiag += diagLocal[Dk];
		}
		Ddiag = A[Dcol + DcolNpad] - Ddiag;
	}
	barrier(CLK_LOCAL_MEM_FENCE);
 
	// diagLocal =  diag * A[Dcol, ]
	// Re-use diagLocal, divide by DL?
	for(Dk = Dlocal; Dk < Dcol; Dk += NlocalSize) {
		DL = A[Dcol + Dk * Npad];
		diagLocal[Dk] = diag[Dk] * DL;
	}
	barrier(CLK_LOCAL_MEM_FENCE); 

	if(N - Dglobal - 1 > Dcol) {
	for(Drow = (N-Dglobal-1); Drow > Dcol; Drow -= Nsize) {

		DL = A[Drow + DcolNpad];

		DrowDk=Drow;
		for(Dk = 0; Dk < Dcol; Dk++) {
//			DL -= A[DrowDk] * diagLocal[Dk];
			DrowDk += Npad;
			DL -= A[Drow + Dk * Npad] * diagLocal[Dk];
			// DL -= A[Drow + Dk * Npad] * A[Dcol + DkNpad] * diag[Dk];
		} // Dk
		A[Drow + DcolNpad] = DL / Ddiag;

	}} // Drow
	// Ddiag is now diag[Dcol]
	// copy it to global memory
	if(Dglobal == 0) {
		diag[Dcol] = Ddiag;
		A[Dcol + DcolNpad] = 1.0;
	}



}


