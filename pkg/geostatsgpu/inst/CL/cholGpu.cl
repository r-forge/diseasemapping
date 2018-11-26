
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

	unsigned int Drow, Dk, DrowDk;

	__local double Ddiag;
	double DL;

	// compute diag[Dcol]
	if(Dlocal == 0) {
		Ddiag = A[Dcol + DcolNpad];
	}
	barrier(CLK_LOCAL_MEM_FENCE);

	for(Dk = Dlocal; Dk < Dcolm1; Dk += NlocalSize) {
		// copy diag to local memory
		diagLocal[Dk] = diag[Dk];
		DL = A[Dk + DcolNpad];
		DL = DL * DL * diagLocal[Dk];
		Ddiag -= DL;
	}
	barrier(CLK_LOCAL_MEM_FENCE);


	// diagLocal =  diag * A[Dcol, ]
	for(Dk=Dlocal; Dk < Dcol; Dk+= NlocalSize) {
		diagLocal[Dk] *= A[Dcol+ Dk * Npad];
	}
	barrier(CLK_LOCAL_MEM_FENCE); 


	for(Drow = N-Dglobal-1; Drow > Dcol; Drow -= Nsize) {

		DL = A[Drow + DcolNpad];

		DrowDk=Drow;
		for(Dk = 0; Dk < Dcolm1; Dk++) {
			DL -= A[DrowDk] * diagLocal[Dk];
			DrowDk += Npad;
			// DL -= A[Drow + Dk * Npad] * diagLocal[Dk];
			// DL -= A[Drow + Dk * Npad] * A[Dcol + DkNpad] * diag[Dk];
		} // Dk
//		A[Drow + DcolNpad] = DL / Ddiag;
	} // Drow
	// Ddiag is now diag[Dcol]
	// copy it to global memory
	if(Dglobal == 0) {
		diag[Dcol] = Ddiag;
	}
if(Dcol == N-1) {
	if(Dglobal < N) {
	A[Dglobal] = Dglobal;
	A[Npad + Dglobal] = Dlocal;
//	A[2*Npad + Dglobal] = Dlocal;
//	A[3*Npad + Dglobal] = Dgroup;
//	A[4*Npad + Dglobal] = Dk;
}
}

}


