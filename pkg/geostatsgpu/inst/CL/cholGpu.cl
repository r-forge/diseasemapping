
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
	const int NinLastCycle,
	const int rowToStartGroupwise,
	const int NlocalSum// typicaly 4
	) {

	const int Nsize = get_global_size(0);
	const int NlocalSize = get_local_size(0);

	const int Dlocal = get_local_id(0);
	const int Dgroup = get_group_id(0);
	const int Ngroup = Nsize / NlocalSize;
	const int DrowStart = Dcol+1+get_global_id(0);

	int Dcycle, DcycleNlocalStorage, Drow, Dk;
	double DL, diagLocalStore;
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
		// do one row per work item
		for(Drow = DrowStart; Drow < rowToStartGroupwise; Drow+= Nsize) {

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
		// loop through remaining rows
		// do one row per work group
		// the first Nlocal elements of diagLocal will be overwritten
		diagLocalStore = diagLocal[Dlocal];
		for(Drow = rowToStartGroupwise + Dgroup; Drow < N; Drow+= Ngroup) {
			// assuming Dlocal < NlocalStorage
			// or more local storage than local work items

			Dk = Dlocal;
			DL = A[Drow + 
					(Dk + DcycleNlocalStorage) * Npad
					] * diagLocalStore;

			for(Dk = Dlocal + NlocalSize; Dk < NlocalStorage;
				Dk += NlocalSize) {
				DL += A[
					Drow + 
					(Dk + DcycleNlocalStorage) * Npad
					] * diagLocal[Dk];
			}
			diagLocal[Dlocal] = DL;
			barrier(CLK_LOCAL_MEM_FENCE);
			// add the local storage
			if(Dlocal < NlocalSum) {
				DL = 0.0;
				for(Dk = Dlocal; Dk < NlocalSize; 
					Dk += NlocalSum) {
					DL += diagLocal[Dk];
				}
				diagLocal[Dlocal] = DL;
			}
			barrier(CLK_LOCAL_MEM_FENCE);

			if(Dlocal==0) {
				DL = 0.0;
				for(Dk = 0; Dk < NlocalSum; Dk++) {
					DL += diagLocal[Dk];
				}
				A[Drow + DcolNpad] -= DL;
			}

		}
	} // end Dcycle
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
	for(Drow = DrowStart; Drow < rowToStartGroupwise; Drow+= Nsize) {

		DL = A[Drow + DcolNpad];
		for(Dk = 0; Dk < NinLastCycle; Dk++) {
			DL -= A[
				Drow + (Dk + NbeforeLastCycle) * Npad
				] * diagLocal[Dk];
		}
		A[Drow + DcolNpad] = DL / diagDcol;
	}

	diagLocalStore = diagLocal[Dlocal];
	for(Drow = rowToStartGroupwise + Dgroup; Drow < N; Drow += Ngroup) {
		if(Dlocal < NinLastCycle) {
			Dk = Dlocal;
			DL = A[Drow + 
					(Dk + DcycleNlocalStorage) * Npad
					] * diagLocalStore;
		} else {
			DL = 0.0;
		}
		for(Dk = Dlocal + NlocalSize; Dk < NinLastCycle;
			Dk += NlocalSize) {
			DL += A[
				Drow + 
				(Dk + DcycleNlocalStorage) * Npad
				] * diagLocal[Dk];
		}
		diagLocal[Dlocal] = DL;
		barrier(CLK_LOCAL_MEM_FENCE);
		// add the local storage
		if(Dlocal < NlocalSum) {
			DL = 0.0;
			for(Dk = Dlocal; Dk < NlocalSize; 
				Dk += NlocalSum) {
				DL += diagLocal[Dk];
			}
			diagLocal[Dlocal] = DL;
		}
		barrier(CLK_LOCAL_MEM_FENCE);

		if(Dlocal==0) {
			DL = 0.0;
			for(Dk = 0; Dk < NlocalSum; Dk++) {
				DL += diagLocal[Dk];
			}
			A[Drow + DcolNpad] -= DL;
		}

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
	const int Dgroup = get_group_id(0);
//	const int Ngroups = Nglobal / Nlocal;
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



		// Drow1 loops through rows after the Dcol'th
		for(Drow1 = 0; Drow1 < Nrows; Drow1++) {
			// index for crossprodLocal[1,Drow1,Dlocal]
			Drow1Nrows = crossprodIndex + Drow1 * Nrows;

			// A[Drow1, Dk] * diag[Dk]
			ADhere = A[DkNpad + Drow1] * Dhere;

			// Drow2 is 'columns' of the crossproduct
			for(Drow2 = Drow1; Drow2 < Nrows; Drow2++){
				// c[Drow1, Drow1, Dlocal] = 
				//   A[Drow1, Dk] * diag[Dk] * A[Drow2,Dk]
				crossprodLocal[Drow1Nrows + Drow2] =
					ADhere * A[DkNpad + Drow2]; 
			}
		}

	} else {
		// Dk >= Dcol, set crossprodLocal to zero
		DkNpad = Dk * Npad + Dcol;

		// Drow1 loops through rows after the Dcol'th
		for(Drow1 = 0; Drow1 < Nrows; Drow1++) {
			// index for crossprodLocal[1,Drow1,Dlocal]
			Drow1Nrows = crossprodIndex + Drow1 * Nrows;
			// Drow2 is 'columns' of the crossproduct
			for(Drow2 = Drow1; Drow2 < Nrows; Drow2++){
				// c[Drow1, Drow1, Dlocal] = 
				//   A[Drow1, Dk] * diag[Dk] * A[Drow2,Dk]
				crossprodLocal[Drow1Nrows + Drow2] = 0.0;
			}
		}
	}

	// now loop through the rest of the columns
	for(Dk = Dglobal+ Nglobal; Dk < Dcol; Dk += Nglobal) {
		Dhere = diag[Dk];
		// index for A[Dcol, Dk]
		DkNpad = Dk * Npad + Dcol;
	
		for(Drow1 = 0; Drow1 < Nrows; Drow1++) {
//			Drow1Nrows = Drow1 * Nrows;
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

	if(Dlocal < NlocalSum ) {

		DkNpad = Dlocal * NrowsSq;
		for(Dk = Dlocal + NlocalSum; Dk < Nlocal; Dk += NlocalSum) {
			// index for 
			// cl[1,1,Dk]
			DlocalNsq = Dk * NrowsSq;
			for(Drow2 = 0; Drow2 < Nrows; Drow2++){
				Drow1Nrows = Drow2 * Nrows;
				for(Drow1 = Drow2; Drow1 < Nrows; Drow1++) {

					crossprodLocal[
						Drow1Nrows + Drow1 + DkNpad //Dlocal * NrowsSq
						] +=
					crossprodLocal[
						Drow1Nrows + Drow1 + DlocalNsq//Dk*NrowsSq
						];

				}
			}
		}
	}

	// sum over the local work items
	if(Dlocal == 0 ) {

#ifdef DEBUG
		for(Dk = 1; Dk < NrowsSq; Dk++) {
		crossprod[
				Dk + Dgroup * NrowsSq
				] = crossprodLocal[Dk];
		}
#endif

		for(Dk = 1; Dk < NlocalSum; Dk ++) {
			// index for 
			// cl[1,1,Dk]
			DlocalNsq = Dk * NrowsSq;


			for(Drow2 = 0; Drow2 < Nrows; Drow2++){
				Drow1Nrows = Drow2 * Nrows;
				for(Drow1 = Drow2; Drow1 < Nrows; Drow1++) {

					crossprodLocal[
						Drow1Nrows + Drow1
						] +=
					crossprodLocal[
						Drow1Nrows+ Drow1 + DlocalNsq//Dk*NrowsSq
						];

				}
			}
		}

		DlocalNsq = Dgroup * NrowsSq;

		for(Drow2 = 0; Drow2 < Nrows; Drow2++){
			Drow1Nrows = Drow2 * Nrows;
			Dk = Drow1Nrows + DlocalNsq;
			for(Drow1 = 0; Drow1 < Nrows; Drow1++) {

					crossprod[
//						Dgroup * NrowsSq +
//						Drow1Nrows
						Dk + Drow1
						] =
					crossprodLocal[
						Drow1Nrows + Drow1
						];

			}
		}
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

		for(Drow2 = Dlocal; //Drow1 + Dlocal; 
			Drow2 < Nrows; 
			Drow2 += Nlocal) {

			theSum = 0.0;
			
			// index for c[Drow2, Drow1, 1]		
			DlIndex = Drow1 * Nrows + Drow2;

			// Dsum is over 3rd dimension
			for(Dsum = 0; Dsum < NtoAdd; Dsum++) {

				// theSum += crossprod[Drow2, Drow1, Dsum]
				theSum += crossprod[
					DlIndex + Dsum * NrowsSq];
//					Drow1*Nrows + Drow2 + Nrows*Nrows*Dsum];
			}
			
			crossprod[DlIndex] = theSum;
//			crossprod[Drow1 * Nrows + Drow2] = theSum;
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

	// crossprod is NrowsCrossprod by NrowsCrossprod
	const int NrowsCrossprod,

	// colSinceCrossprod * (NrowsCrossprod + 1)
	const int colSinceCrossprodNrowsCrossprodP1, 

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
		Dk < NrowsCrossprod;
		Dk+=Nlocal) {

		diagLocal[Dk] = A[
//			(Dk + colAtCrossprod )*Npad + Dcol
				AforCrossprodStartJ + Dk * Npad
				] * diag[colAtCrossprod + Dk];
	}
	barrier(CLK_LOCAL_MEM_FENCE);

	// update the A[,Dcol] to
	// A[Drow, Dcol] = sum_k L_ik L_jk D_k

	for(Drow = Dglobal; 
		Drow < Nrows; 
		Drow += Nglobal) {

		Dcrossprod = crossprod[
//			colSinceCrossprod * (NrowsCrossprod + 1) +
			colSinceCrossprodNrowsCrossprodP1 + 
			Drow];

#ifdef DEBUG
		A[colSinceCrossprod * Npad + 200 + Drow] = Dcrossprod;
		Dk = 0;
		A[colSinceCrossprod * Npad + 400 + Drow] = A[
					(Dk + colAtCrossprod )*Npad + 
					Dcol];
		A[colSinceCrossprod * Npad + 600 + Drow] = 				A[
					(Dk + colAtCrossprod )*Npad + 
					Dcol + Drow + Dk
				];
		A[colSinceCrossprod * Npad + 800 + Drow] = 				diag[colAtCrossprod + Dk];
#endif


		AforCrossprodStartI = 
			AforCrossprodStartJ + Drow;

		for(Dk = 0; Dk < colSinceCrossprod; Dk++) {
			Dcrossprod +=  
				A[
//					(Dk + colAtCrossprod )*Npad + 
//					Dcol + Drow
				AforCrossprodStartI + Dk * Npad
				]  * diagLocal[Dk];
//				A[
//					(Dk + colAtCrossprod )*Npad + 
//					Dcol] * 
//				diag[colAtCrossprod + Dk];
		}

#ifdef DEBUG
		A[colSinceCrossprod * Npad + 1000 + Drow] = Dcrossprod;
#endif




// index for A[DrowFromZero, Dcol]
//  DrowFromZero = Drow + Dcol
//		Dk = AforUpdateStart + Drow;

		if(Drow == 0) {
			diag[Dcol] = A[
				AforUpdateStart
//				Dcol*Npad + Dcol
				] - 
			Dcrossprod;
		} else {
			A[
				//Dcol*Npad + Dcol 
				AforUpdateStart + Drow
			] -= Dcrossprod;
		}

		// note havent yet divided by diagonal

	}

}