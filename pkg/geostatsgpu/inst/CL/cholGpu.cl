
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

		for(Dk = Dlocal; Dk < NlocalStorage; 
			Dk += NlocalSize) {
	
			diagLocal[Dk] = diagTimesRowOfA[
				Dk + DcycleNlocalStorage];

		}

		// synchronize all the work items
		// so that all of diagLocal is ready
		barrier(CLK_LOCAL_MEM_FENCE);

		// loop through rows
		for(Drow = DrowStart; Drow < N; Drow+= Nsize) {

			// do this cycle's bit of the cholesky
			DL = A[Drow + DcolNpad];
			for(Dk = 0; Dk < NlocalStorage; Dk++) {
				DL -= 
				A[
					Drow + 
					(Dk + DcycleNlocalStorage) * Npad
					] * 
				diagLocal[Dk];
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
				Dk + NbeforeLastCycle];
				//DcycleNlocalStorage];
	
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
	__global double *result,
	__local double *diagLocal,
	const int N,
	const int NlocalSum// typicaly 4, 
	) 
{


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
	if(Dlocal < NlocalSum) {
		Ddiag = 0.0;
		for(Dk = Dlocal; Dk < NlocalSize; 
			Dk += NlocalSum) {
			Ddiag += diagLocal[Dk];
		}
		diagLocal[Dlocal] = Ddiag;
	}
	barrier(CLK_LOCAL_MEM_FENCE);

	if(Dlocal==0) {
		Ddiag = 0.0;
		for(Dk = 0; Dk < NlocalSum; Dk++) {
			Ddiag += diagLocal[Dk];
		}

		result[get_group_id(0)] = Ddiag;
	}


}

// crossprod_ij contains sum_k L_kj L_ij D_k
// crossprodLocal must contain over Nlocal * Nrows^2 
// entries

__kernel void cholCrossprod(
	__global double *A, 
	__global double *diag,
	// Nrows by Nrows by Ngroups
	__global double *crossprod,
	// Nrows by Nrows by Nlocal
	__local double *crossprodLocal,
	const int Dcol,
	const int DcolNpad,
	const int Nrows, // N - Dcol
	const int NrowsSq,
	const int Npad,
	const int NlocalSum// typicaly 4, 
	) {


	const int Dlocal = get_local_id(0);
	const int Nlocal = get_local_size(0);
	const int Dglobal = get_global_id(0);
	const int Nglobal = get_global_size(0);
// location in crossprodLocal to store the results
	// index for crossprodLocal[1,1,Dlocal]
	const int crossprodIndex = Dlocal*NrowsSq;

	int Drow1, Drow2, Dk, Drow1Nrows; 
	int DkNpad, DlocalNsq, DlIndex, DkIndex;
	double ADhere, Dhere;

	// Dk loops through columns of A

	// do the first column, 
	// to initialize crossprodLocal
	Dk = Dglobal;

	if(Dk < Dcol) {
	Dhere = diag[Dk];

	// index for A[Dcol, Dk]
	DkNpad = Dk * Npad + Dcol;
	// Drow1 loops through rows from Dcol
	for(Drow1 = 0; Drow1 < Nrows; Drow1++) {
		// index for crossprodLocal[1,Drow1,Dlocal]
		Drow1Nrows = crossprodIndex + Drow1 * Nrows;

		// A[Drow1, Dk] * diag[Dk]
		ADhere = A[DkNpad + Drow1] * Dhere;

		for(Drow2 = Drow1; Drow2 < Nrows; Drow2++){
			// c[Drow1, Drow1, Dlocal] = 
			//   A[Drow1, Dk] * diag[Dk] * A[Drow2,Dk]
			crossprodLocal[Drow1Nrows + Drow2] =
				ADhere * A[DkNpad + Drow2]; 
		}
	}

	}


	// now loop through the rest of the columns
	for(Dk = Dglobal+ Nglobal; Dk < Dcol; Dk += Nglobal) {
		Dhere = diag[Dk];
		// index for A[Dcol, Dk]
		DkNpad = Dk * Npad + Dcol;
	
		for(Drow1 = 0; Drow1 < Nrows; Drow1++) {
			Drow1Nrows = Drow1 * Nrows;
			// index for crossprodLocal[1,Drow1,Dlocal]
			Drow1Nrows = crossprodIndex + Drow1 * Nrows;

			// A[Drow1, Dk] * diag[Dk]
			ADhere = A[DkNpad + Drow1] * Dhere;

			for(Drow2 = Drow1; Drow2 < Nrows; Drow2++){
			// c[Drow1, Drow1, Dlocal] = 
			//   A[Drow1, Dk] * diag[Dk] * A[Drow2,Dk]
				crossprodLocal[Drow1Nrows + Drow2] +=
					ADhere * A[DkNpad + Drow2]; 
			}
		}
	}
	barrier(CLK_LOCAL_MEM_FENCE);

	// sum over the local work items
	if(Dlocal < NlocalSum) {
		// index for 
		// crossprodLocal[1,1,Dlocal]
		DlocalNsq = Dlocal * NrowsSq;

		for(Dk = Dlocal + NlocalSum; Dk < Nlocal; 
			Dk += NlocalSum) {
			// index for 
			// crossprodLocal[1,1,Dk]
			DkNpad = Dk * NrowsSq;

			// Drow1 = dim 2 of crossprodLocal
			for(Drow1 = 0; Drow1 < Nrows; Drow1++) {

				Drow1Nrows = Drow1 * Nrows;

				// index for 
				// crossprodLocal[1,Drow1,Dlocal]
				DlIndex = DlocalNsq + Drow1Nrows;	
				// index for 
				// crossprodLocal[1,Drow1,Dk]
				DkIndex = DkNpad + Drow1Nrows;

				for(Drow2 = Drow1; Drow2 < Nrows; Drow2++){

					// cl[D2,D1,Dlocal] +=
					// cl[D2,D1,Dk]

					crossprodLocal[
						DlIndex + Drow2
					] += crossprodLocal[
						DkIndex + Drow2
					];
				}
			}
		}
	}
	barrier(CLK_LOCAL_MEM_FENCE);

	// sum over the local work items
	if(Dlocal == 0 ) {
		for(Dk = 1; Dk < NlocalSum; Dk ++) {
			// index for 
			// cl[1,1,Dk]
			DkNpad = Dk * NrowsSq;

			for(Drow1 = 0; Drow1 < Nrows; Drow1++) {

				Drow1Nrows = Drow1 * NrowsSq;

				// index for 
				// crossprodLocal[1,Drow1,1]
				DlIndex = Drow1Nrows;	
				// index for 
				// crossprodLocal[1,Drow1,Dk]
				DkIndex = DkNpad + Drow1Nrows;

				for(Drow2 = Drow1; Drow2 < Nrows; Drow2++){
					crossprodLocal[
						DlIndex + Drow2
					] += crossprodLocal[
						DkIndex + 
						Drow2
					];
				}
			}
		}
		// copy to global memory
		// address for c[1,1,Dgroup]
		DkNpad = get_group_id(0) * NrowsSq;
	
		for(Drow1 = 0; Drow1 < Nrows; Drow1++) {
			Drow1Nrows = Drow1 * Nrows;	

			// index for c[1,Drow1,Dgroup]
			DlIndex = DkNpad + Drow1Nrows;

			// index for cl[1,Drow1,Dlocal]
			DkIndex = Drow1Nrows;

			for(Drow2 = Drow1; Drow2 < Nrows; Drow2++){
				crossprod[
					DlIndex + Drow2
				] = crossprodLocal[
					DkIndex + Drow2
				];
			}
		}  // end Drow1
	} // end Dlocal == 0

}

// crossprod is Nrows by Nrows by NtoAdd
// apply(crossprod,c(1,2), sum)
__kernel void cholSumCrossprod(
	__global double *crossprod,
	const int Nrows,
	const int NrowsSq,
	const int NtoAdd // number of groups in cholCrossprod
	) {



	const int Dlocal = get_local_id(0);
	const int Dgroup = get_group_id(0);
	const int Ngroups = get_num_groups(0);
	const int Nlocal = get_local_size(0);

	int Drow1, Drow2, Dsum;
	int DlIndex;
	double theSum;
	// groups are dimension 2
	// local items are dimension 1
	for(Drow1 = Dgroup; Drow1 < Nrows; Drow1 += Ngroups) {

		for(Drow2 = Drow1 + Dlocal; Drow2 < Nrows; 
			Drow2 += Nlocal) {

			theSum = 0.0;
			
			// index for c[Drow2, Drow1, 1]		
			DlIndex = Drow1 * Nrows + Drow2;

			// Dsum is over 3rd dimension
			for(Dsum = 0; Dsum < NtoAdd; Dsum++) {

				// theSum += crossprod[Drow2, Drow1, Dsum]
				theSum += crossprod[
					DlIndex + Dsum * NrowsSq];
			}
			
			crossprod[DlIndex] = theSum;
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

	// index for A[Dcol, Dcol]
	const int AforUpdateStart,

	// crossprod is NrowsC by NrowsC
	const int NrowsCrossprod, // N - colAtCrossprod

	// number of rows of A which need updating
	const int Nrows, // N - Dcol

	const int Npad
	) {


	const int Nglobal = get_global_size(0);
	const int Dglobal = get_global_id(0);
	const int Nlocal = get_local_size(0);

	double Dcrossprod;
	int Drow, Dk, AforCrossprodStartI;

	// compute the L_jk D_k
	// where k = colAtCrossprod
	// store locally
	for(Dk = get_local_id(0); 
		Dk < colSinceCrossprod;
		Dk+=Nlocal) {

		diagLocal[Dk] = A[
				AforCrossprodStartJ + Dk * Npad
				] * diag[colAtCrossprod + Dk];
	}
	barrier(CLK_LOCAL_MEM_FENCE);


	// update the A[,Dcol] to
	// A[Drow, Dcol] = sum_k L_ik L_jk D_k

	for(Drow = Dglobal; Drow < Nrows; 
		Drow += Nglobal) {

		// index for A[Drow, Dcol]
		AforCrossprodStartI = 
			AforCrossprodStartJ + Drow;

		// stored crossprod = 
		// sum_{k=1}^colsSinceCrossprod
		//    A[Drow,k] * A[Dcol,k] * D_k	
		Dcrossprod = crossprod[
			colSinceCrossprod + Drow];

		// update with
		// sum_{k=colsSince + 1}^{j-1}
		//    A[i,k] * A[j,k] * D_k	
		for(Dk = 0; Dk < colSinceCrossprod; Dk++) {
			Dcrossprod += A[
				AforCrossprodStartI + Dk * Npad
				] * diagLocal[Dk];
		}

		// index for A[DrowFromZero, Dcol]
		//  DrowFromZero = Drow + Dcol
		Dk = AforUpdateStart + Drow;
		if(Drow == 0) {
			diag[Dcol] = (A[Dk] - Dcrossprod);
		} else {
			A[Dk] = (A[Dk] - Dcrossprod);
		}
		// note havent yet divided by diagonal

	}

}