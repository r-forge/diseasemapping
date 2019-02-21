
// https://en.wikipedia.org/wiki/Cholesky_decomposition#LDL_decomposition_2


// crossprod_ij contains sum_k L_kj L_ij D_k
// entries where k=Diter are new

__kernel void cholUpdateCrossprod(
	__global double *A, 
	__global double *crossprod,
	const double Ddiag,
	const int Diter,
	const int DiterNpad,
	const int N, const int Npad
	) {

	int Drow, Drow, Dcell;
	const int Dlocal = get_local_id(0);
	const int Nlocal = get_local_size(0);
	const int Dgroup = get_group_id(0);
	const int Ngroup = get_group_size(0);
	__local double ADrowDiter;
	// each row a work group
	// each column a work item

	for(Drow = Diter + Dgroup; Drow < N; 
		Drow += Ngroup) {
		if(Dlocal == 0) {
			ADrowDiter = A[Drow + DiterNpad];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	
		for(Dcol = 0; Dcol < Diter; ++Dcol) {
			Dcell = Dcol*Npad + Drow;

			crossprod[Dcell] += ADrowDiter *
				A[Dcol + DiterNpad] * Ddiag;
		}

	}

}

__kernel void cholUpdateA(){}
