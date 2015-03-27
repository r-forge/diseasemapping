#include<R.h>
#include<Rmath.h>
#include<R_ext/Lapack.h>
#include<R_ext/Applic.h>
#include<R_ext/Print.h>
#include<R_ext/Utils.h>
//#include"/usr/lib/R/library/Matrix/include/Matrix.h"
#include<R_ext/Rdynload.h>
#include<Matrix.h>
#include<Matrix_stubs.c>



int cholmod_solve2         /* returns TRUE on success, FALSE on failure */
(
    /* ---- input ---- */
    int sys,		            /* system to solve */
    cholmod_factor *L,	            /* factorization to use */
    cholmod_dense *B,               /* right-hand-side */
    cholmod_sparse *Bset,
    /* ---- output --- */
    cholmod_dense **X_Handle,       /* solution, allocated if need be */
    cholmod_sparse **Xset_Handle,
    /* ---- workspace  */
    cholmod_dense **Y_Handle,       /* workspace, or NULL */
    cholmod_dense **E_Handle,       /* workspace, or NULL */
    /* --------------- */
    cholmod_common *Common
);

int attribute_hidden Mbob_cholmod_solve2(
		int sys,
		CHM_FR L,
		CHM_DN B,//right
		CHM_DN *X,//solution
		CHM_DN *Yworkspace,
		CHM_DN *Eworkspace,
		cholmod_common *c)
{
    static int(*fun)(
    		int,
    		const_CHM_FR, // L
    		const_CHM_DN, // B
    		CHM_SP, // Bset
			CHM_DN*, // X
			CHM_DN*, // Xset
			CHM_DN*, // Y
			CHM_DN*, // E
			cholmod_common*) = NULL;

    if (fun == NULL)
    	fun = (int(*)(int,
    		const_CHM_FR, // L
    		const_CHM_DN, // B
    		CHM_SP, // Bset
			CHM_DN*, // X
			CHM_DN*, // Xset
			CHM_DN*, // Y
			CHM_DN*, // E
			cholmod_common*)
	    )R_GetCCallable("Matrix", "cholmod_solve2");

    return fun(
    		sys,
    		L,
			B,NULL,
			X, NULL,
			Yworkspace,
			Eworkspace,
			c);
}

// global variables, for when this is in an optimizer
CHM_SP Q;
CHM_FR L;
CHM_DN obsCovRot, Lx, DLx;
CHM_DN YwkL, EwkL, YwkD, EwkD; // workspaces
cholmod_common c;
double logLtwo[2], detTwo[2], *YXVYXglobal, *YXYX, detQ;
int Nxy, Nobs, Ncov, Nxysq;
int Ltype;

// compute sums of squares from cross products
void ssqFromXprod(
		double *YXVinvYX, // N by N
		double *detXVinvX,
		int N
){

	int oneI=1, infoCholXX, infoInvXX, Np1, Nm1;
	int D;
	double	oneD=1.0, zeroD=0.0;
	double xybeta;

	Np1 = N+1;
	Nm1 = N-1;

	// invert X Vinv X
	// first cholesky

	//  cholesky X Vinv X
	F77_CALL(dpotrf)("L",
		&Nm1, &YXVinvYX[Np1],
		&N, // NcolM1 by NcolM1 submatrix of Ncol by Ncol matrix
		&infoCholXX);

	*detXVinvX  = 0.0;
	for(D=1;D<N;++D){
		*detXVinvX  += log(YXVinvYX[D*N+D]);
	}
// then invert
F77_NAME(dpotri)("L",
		&Nm1,
		&YXVinvYX[Np1],
		&N,
		&infoInvXX);

// put beta hat in first column (first row still has LyLx)
F77_NAME(dsymm)(
		"R", "L", // A on right, A in lower
		&oneI, &Nm1, // C has one row, Ncov columns
		&oneD,
		&YXVinvYX[Np1], &N,// XVinvX^(-1)
		&YXVinvYX[N],&N, // Lx
		&zeroD,
		&YXVinvYX[1], &oneI //betahat
		);

// LxLy beta
F77_NAME(dgemm)(
		//	op(A), op(B),
		"N", "N",
		// nrows of op(A), ncol ob(B), ncol op(A) = nrow(opB)
  		&oneI, &oneI, &Nm1,
		// alpha
		&oneD,
		// A, lda(A)
		&YXVinvYX[Np1], &N,
		// B, lda(B)
		&YXVinvYX[1], &Nm1,
		// beta
  		&zeroD,
		// C, nrow(c)
		&xybeta, &oneI);

// total ssq
YXVinvYX[0] -= xybeta;

}

// logL given xisqTausq
// needs global variables
// Q, L, c, detTwo
// obsCovRot, Lx, YwkL, EwkL, DLx, YwkD, EwkD
// YXVYXglobal, YXYX, Nxy, Nobs,
double logLoneNugget(double xisqTausq){

	double minusXisqTausq, zeroD=0.0, moneD=-1.0;
	double *DYXVYX, result;
	int oneI=1;

	DYXVYX=YXVYXglobal;

	M_cholmod_factorize_p(
		Q,
		&xisqTausq, // beta
		(int*)NULL, 0 /*fsize*/,
		L, &c
	);

	detTwo[0] = M_chm_factor_ldetL2(L);
//Lx =
cholmod_solve2(
		CHOLMOD_L,
		L,
		obsCovRot, NULL,
		&Lx, NULL,
		&YwkL, &EwkL,
		&c);
//DLx =
cholmod_solve2(
		CHOLMOD_D,
		L,
		Lx, NULL,
		&DLx, NULL,
		&YwkD, &EwkD,
		&c);

// cross product
minusXisqTausq = -xisqTausq;

// C := alpha*op( A )*op( B ) + beta*C,
F77_NAME(dgemm)(
		//	op(A), op(B),
		"T", "N",
		// nrows of op(A), ncol ob(B), ncol op(A) = nrow(opB)
  		&Nxy, &Nxy, &Nobs,
		// alpha
		&minusXisqTausq,
		// A, nrow(A)
  		Lx->x, &Nobs,
		// B, nrow(B)
		DLx->x, &Nobs,
		// beta
  		&zeroD,
		// C, nrow(c)
		DYXVYX, &Nxy);

// add the cross prod of data
F77_NAME(daxpy)(
		&Nxysq,
		&moneD,
		YXYX,
		&oneI,
		DYXVYX,
		&oneI);

ssqFromXprod(
		DYXVYX, // N by N
		&detTwo[1],
		Nxy);

// ml
logLtwo[0] = Nobs*log(DYXVYX[0]/Nobs) - detQ + 2*detTwo[0];
// reml
logLtwo[1] = (Nobs-Nxy+1)*log(DYXVYX[0]/(Nobs-Nxy+1)) -
		detQ + 2*detTwo[0] - 2*detTwo[1];

if(Ltype){
	result = logLtwo[1];
} else {
	result = logLtwo[0];
}
return result;

}

SEXP gmrfLik(
		SEXP QR,
		SEXP obsCovR,
		SEXP xisqTausq
		){

	int DxisqTausq, NxisqTausq; // length(xisqTausq)
	double	oneD=1.0, zeroD=0.0;
	double *YXVYX, *determinant, *determinantForReml;
	double *m2logL, *m2logReL;
	SEXP resultR;
	CHM_DN obsCov;

	Ltype=0; // set to 1 for reml

	Nobs = INTEGER(GET_DIM(obsCovR))[0];
	Nxy = INTEGER(GET_DIM(obsCovR))[1];
	NxisqTausq = LENGTH(xisqTausq);
	Nxysq = Nxy*Nxy;



	resultR = PROTECT(allocVector(REALSXP, (4+Nxysq)*NxisqTausq));
	YXVYX = REAL(resultR);
	YXYX = (double *) calloc(Nxysq,sizeof(double));
	determinant = &YXVYX[Nxysq*NxisqTausq];
	determinantForReml = &YXVYX[(1+Nxysq)*NxisqTausq];
	m2logL = &YXVYX[(2+Nxysq)*NxisqTausq];
	m2logReL = &YXVYX[(3+Nxysq)*NxisqTausq];

	Q = AS_CHM_SP(QR);
	obsCov = AS_CHM_DN(obsCovR);


	M_R_cholmod_start(&c);

	// get some stuff ready
	// Q = P' L D L' P
	L = M_cholmod_analyze(Q, &c);
	M_cholmod_factorize(Q,L, &c);

	obsCovRot = M_cholmod_solve(CHOLMOD_Pt, L,obsCov,&c);




	// likelihood without nugget

	// determinant
	determinant[0] = M_chm_factor_ldetL2(L);
	detQ = determinant[0];


	// cross product of data
	F77_NAME(dgemm)(
			//	op(A), op(B),
			"T", "N",
			// nrows of op(A), ncol ob(B), ncol op(A) = nrow(opB)
	  		&Nxy, &Nxy, &Nobs,
			// alpha
			&oneD,
			// A, nrow(A)
			obsCovRot->x, &Nobs,
			// B, nrow(B)
			obsCovRot->x, &Nobs,
			// beta
	  		&zeroD,
			// C, nrow(&c)
			YXYX, &Nxy);


	// allocate Lx
	Lx = M_cholmod_copy_dense(obsCovRot,&c);

	// YX Vinv YX
		M_cholmod_sdmult(
				Q,
				0, &oneD, &zeroD, // transpose, scale, scale
				obsCov,Lx,// in, out
				&c);




	// put t(obscov) Q obscov in result
	F77_NAME(dgemm)(
			//	op(A), op(B),
			"T", "N",
			// nrows of op(A), ncol ob(B), ncol op(A) = nrow(opB)
	  		&Nxy, &Nxy, &Nobs,
			// alpha
			&oneD,
			// A, nrow(A)
	  		Lx->x, &Nobs,
			// B, nrow(B)
			obsCov->x, &Nobs,
			// beta
	  		&zeroD,
			// C, nrow(c)
			YXVYX, &Nxy);

	ssqFromXprod(
			YXVYX, // N by N
			determinantForReml,
			Nxy);

	m2logL[0] = Nobs*log(YXVYX[0]/Nobs) - determinant[0];

	m2logReL[0] = (Nobs-Ncov)*log(YXVYX[0]/(Nobs-Ncov)) +
			determinantForReml[0]- determinant[0];


	// now with xisqTausq
	for(DxisqTausq=1;DxisqTausq < NxisqTausq;++DxisqTausq){

		YXVYXglobal = &YXVYX[DxisqTausq*Nxysq];

		logLoneNugget(REAL(xisqTausq)[DxisqTausq]);

		// assign global values into their correct spot
		determinant[DxisqTausq]=detTwo[0];
		determinantForReml[DxisqTausq]=detTwo[1];
		m2logL[DxisqTausq] = logLtwo[0];
		m2logReL[DxisqTausq] = logLtwo[1];

	}

	M_cholmod_free_factor(&L, &c);
	M_cholmod_free_dense(&obsCovRot, &c);

	M_cholmod_free_dense(&Lx, &c);
	M_cholmod_free_dense(&DLx, &c);

// don't free Q because it's from an R object
//	M_cholmod_free_sparse(&Q, &c);

// don't free obsCov because it's from an R object
//	M_cholmod_free_dense(&obsCov, &c);

	free(YXYX);
	M_cholmod_free_dense(&YwkL, &c);
	M_cholmod_free_dense(&YwkD, &c);
	M_cholmod_free_dense(&EwkL, &c);
	M_cholmod_free_dense(&EwkD, &c);

	M_cholmod_finish(&c);

	UNPROTECT(1);
	return resultR;
}


