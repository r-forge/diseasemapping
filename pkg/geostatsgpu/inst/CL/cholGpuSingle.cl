__kernel void cholGpu(
	__global double *A, 
	__global double *diag,
	__local double *diagLocal,
	const int Dcol,
	const int N, 
	const int Npad,
	const int NlocalStorage)
{

	const int Dlocal = get_local_id(0);
	const int Dglobal = get_global_id(0);
	const int Nsize = get_global_size(0);
	const int NlocalSize = get_local_size(0);
 
	const int Dcolm1 = Dcol - 1;
	const int DcolNpad = Dcol * Npad;
	const int maxDcolNlocalStorage = min(Dcol, NlocalStorage);

	__private int Drow, Dk, DrowDk;
	__private double DL, Ddiag;


	Ddiag = 0.0;
	for(Dk = Dglobal; Dk < Dcol; Dk += Nsize) {
		DL = A[Dcol + Dk * Npad];
		DL *= DL;
		diagLocal[Dlocal] += diag[Dk] * DL;
	}
	barrier(CLK_LOCAL_MEM_FENCE);
// add 'em up to get diagonal
	if(Dlocal==0) {
		Ddiag = 0.0;
		for(Dk = 0; Dk < NlocalSize; Dk++) {
			Ddiag += diagLocal[Dk];
		}
		Ddiag = A[Dcol + DcolNpad] - Ddiag;
	}
	barrier(CLK_LOCAL_MEM_FENCE);
 
	// diagLocal =  diag * A[Dcol, ]
	// TO DO: use less memory for diagLocal
	for(Dk = Dlocal; Dk < maxDcolNlocalStorage; Dk += NlocalSize) {
		DL = A[Dcol + Dk * Npad];
		diagLocal[Dk] = diag[Dk] * DL;
	}
	barrier(CLK_LOCAL_MEM_FENCE); 

	for(Drow = (N-Dglobal-1); Drow > Dcol; Drow -= Nsize) {

		DL = A[Drow + DcolNpad];

		DrowDk=Dcolm1*Npad;
		for(Dk = Dcolm1; Dk >= 0; Dk--) {
//			DL -= A[DrowDk] * diagLocal[Dk];
			DrowDk -= Npad;
			DL -= A[Drow + Dk * Npad] * diagLocal[Dk];
			// DL -= A[Drow + Dk * Npad] * A[Dcol + DkNpad] * diag[Dk];
		} // Dk
		A[Drow + DcolNpad] = DL / Ddiag;

	} // Drow
	// Ddiag is now diag[Dcol]
	// copy it to global memory
	if(Dglobal == 0) {
		diag[Dcol] = Ddiag;
		A[Dcol + DcolNpad] = 1.0;
	}


}


