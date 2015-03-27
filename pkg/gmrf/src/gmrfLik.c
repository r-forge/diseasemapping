#include<R.h>
#include<Rmath.h>
#include<R_ext/Lapack.h>
#include<R_ext/Applic.h>
#include<R_ext/Print.h>
#include<R_ext/Utils.h>
#include<Matrix.h>
#include<R_ext/Rdynload.h>

int attribute_hidden
M_cholmod_solve2(
		int sys,
		CHM_FR L,
		CHM_DN B,//right
		CHM_DN *X,//solution
		CHM_DN *Yworkspace,
		CHM_DN *Eworkspace,
		CHM_CM Common)
{
    static int(*fun)(
    		int,
    		const_CHM_FR, // L
    		const_CHM_DN, // B
    		CHM_SP, // Bset
			*CHM_DN, // X
			*CHM_DN, // Xset
			*CHM_DN, // Y
			*CHM_DN, // E
			CHM_CM) = NULL;

    if (fun == NULL)
    	fun = (int(*)(int,
    		const_CHM_FR, // L
    		const_CHM_DN, // B
    		CHM_SP, // Bset
			*CHM_DN, // X
			*CHM_DN, // Xset
			*CHM_DN, // Y
			*CHM_DN, // E
			CHM_CM)
	    )R_GetCCallable("Matrix", "cholmod_solve2");

    return fun(
    		sys,
    		L,
			B,NULL,
			X, NULL,
			Yworkspace,
			Eworkspace,
			Common);
}

// global variables, for when this is in an optimizer
CHM_SP Q;
CHM_FR L;
CHM_DN obsCovRot, Lx;
CHM_DN YwkL, EwkL, YwkD, EwkD; // workspaces
CHM_CM Common;
double logLtwo[2], detTwo[2], *YXVYXglobal, YXYX, detQ;
int Nxy, Nobs, Ncov;
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
	for(D=1;D<Ncol;++D){
		*detXVinvX  += log(YXVinvYX[D*N+D]);
	}
// then invert
F77_NAME(dpotri)("L",
		&Nm1,
		&DYXVYX[Np1],
		&N,
		&infoInvXX);

// put beta hat in first column (first row still has LyLx)
F77_NAME(dssymm)(
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
		&xybeta, oneI);

// total ssq
YXVinvYX[0] -= xybeta;

}

// logL given xisqTausq
// needs global variables
// Q, L, Common, detTwo
// obsCovRot, Lx, YwkL, EwkL, DLx, YwkD, EwkD
// YXVYXglobal, YXYX, Nxy, Nobs,
double logLoneNugget(double xisqTausq){

	double minusXisqTausq, zeroD=0.0, moneD=-1.0;
	double *DYXVYX, result;
	int oneI=1;

	DYXVYX=YXVYXglobal;

	M_cholmod_factorize_p(
		Q,
		xisqTausq, // beta
		(int*)NULL, 0 /*fsize*/,
		L, Common
	);

	detTwo[0] = M_chm_factor_ldetL2(L);
//Lx =
M_cholmod_solve2(
		CHOLMOD_L,
		L,
		obsCovRot,
		&Lx,
		&YwkL, &EwkL,
		Common);
//DLx =
M_cholmod_solve2(
		CHOLMOD_D,
		L,
		Lx,
		&DLx,
		&YwkD, &EwkD,
		Common);

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
	double *YXVYX, determinant, determinantForReml;
	double *m2logL, m2logReL;
	int Nxysq;
	SEXP resultR;
	CHM_DN obsCov;
	CHM_SP Lmat;

	Ltype=0; // set to 1 for reml

	Nobs = INTEGER(GET_DIM(obsCovR))[0];
	Nxy = INTEGER(GET_DIM(obsCovR))[1];
	NxiqTausq = LENGTH(xisqTausq);
	Nxysq = Nxy*Nxy;
	resultR = PROTECT(allocVector(REALSXP), (4+Nxysq)*NxisqTausq);
	YXVYX = REAL(resultR);
	YXYX = (double *) calloc(Nxysq,sizeof(double));
	determinant = &YXVYX[Nxysq*NxisqTausq];
	determinantForReml = &YXVYX[(1+Nxysq)*NxisqTausq];
	m2logL = &YXVYX[(2+Nxysq)*NxisqTausq];
	m2logReL = &YXVYX[(3+Nxysq)*NxisqTausq];

	M_R_cholmod_start(Common);

	Q = AS_CHM_SP(QR);
	obsCov = AS_CHM_DN(obsCovR);

	// get some stuff ready
	// Q = P' L D L' P
	L = M_R_chomod_analyze(Q, Common);
	M_R_cholmod_factorize(Q,L, Common);

	obsCovRot = M_R_cholmod_solve(CHOLMOD_Pt, L,obsCov,Common);

	// cross product of data
	F77_NAME(dgemm)(
			//	op(A), op(B),
			"T", "N",
			// nrows of op(A), ncol ob(B), ncol op(A) = nrow(opB)
	  		&Nxy, &Nxy, &Nobs,
			// alpha
			&oneD,
			// A, nrow(A)
			obsCovRot, Nobs,
			// B, nrow(B)
			obsCovRot, Nobs,
			// beta
	  		&zeroD,
			// C, nrow(c)
			YXYX, Nxy);

	// likelihood without nugget

	// determinant
	determinant[0] = M_chm_factor_ldetL2(L);
	detQ = determinant[0];

	Lmat = M_cholmod_factor_to_sparse(L, Common);

	M_cholmod_sdmult(
			Lmat,
			1, 1.0, 0.0, // transpose, scale, scale
			obsCovRot,Lx,// in, out
			Common);
	DLx = M_R_cholmod_solve(CHOLMOD_D, L, Lx, Common);

	// put t(obscov) Q obscov in result
	F77_NAME(dgemm)(
			//	op(A), op(B),
			"T", "N",
			// nrows of op(A), ncol ob(B), ncol op(A) = nrow(opB)
	  		&Nxy, &Nxy, &Nobs,
			// alpha
			&oneD,
			// A, nrow(A)
	  		Lx->x, Nobs,
			// B, nrow(B)
			DLx->x, Nobs,
			// beta
	  		&zeroD,
			// C, nrow(c)
			YXVYX, Nxy);

	ssqFromXprod(
			YXVYX, // N by N
			detForReml,
			Nxy);

	m2logL[0] = Nobs*log(YXVYX[0]/Nobs) - determinant[0];

	m2logReL[0] = (Nobs-Ncov)*log(YXVYX[0]/(Nobs-Ncov)) +
			detForReml[0]- determinant[0];

	// now with xisqTausq
	for(DxisqTausq=1;DxisqTausq < NxiqTausq;++DxisqTausq){

		YXVYXglobal = &YXVYX[DxisqTausq*Nxysq];

		logLoneNugget(REAL(xisqTausq)[DxisqTausq]);

		// assign global values into their correct spot
		determinant[DxisqTausq]=detTwo[0];
		detForReml[DxisqTausq]=detTwo[1];
		m2logL[DxisqTausq] = logLtwo[0];
		m2logReL[DxisqTausq] = logLtwo[1];

	}


	M_cholmod_free_dense(Lx, Common);
	M_cholmod_free_dense(DLx, Common);
	free(YXYX);
	M_cholmod_free_dense(YwkL, Common);
	M_cholmod_free_dense(YwkD, Common);
	M_cholmod_free_dense(EwkL, Common);
	M_cholmod_free_dense(EwkD, Common);
	M_cholmod_free_dense(obsCov, Common);
	M_cholmod_free_dense(obsCovRot, Common);
	M_cholmod_free_sparse(Q, Common);
	M_cholmod_free_factor(L, Common);
	M_cholmod_free_sparse(Lmat, Common);
	M_cholmod_finish(&Common);

	UNPROTECT(1);
	return resultR;
}


