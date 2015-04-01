/*
 * needs modified matrix package
 * svn checkout svn://svn.r-forge.r-project.org/svnroot/matrix/pkg/Matrix_1.1-5-branch
 * add line     RREGDEF(cholmod_solve2);
 *  to  Matrix_1.1-5-branch/src/init.c
 */

#include<R.h>
#include<Rmath.h>
#include<R_ext/Lapack.h>
#include<R_ext/Applic.h>
#include<R_ext/Print.h>
#include<R_ext/Utils.h>
#include<R_ext/Rdynload.h>
#include<Matrix.h>
#include<Matrix_stubs.c>

/*
	a function for interfacing to the Matrix package
	hopefully this will be in Matrix_stubs.c soon
		*/
attribute_hidden int M_cholmod_solve2(
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

/*
 *  global variables, so that an optimizer can be build eventually
 */
CHM_SP Q;
CHM_FR L;
CHM_DN obsCovRot, Lx;
CHM_DN YwkL, EwkL, YwkD, EwkD; // workspaces
cholmod_common c;
double *logLtwo, detTwo[2], *YXVYXglobal, *YXYX, *YrepAdd;
double *copyLx;
int Nxy, Nobs, Ncov, Nrep, Nxysq;
int Ltype;

/*
 * compute sums of squares from cross products
 *
 */
void ssqFromXprod(
		double *YXVinvYX, // N by N
		double *detXVinvX,
		const int N, const int Nrep,
		double *copyLx
){

	int oneI=1, infoCholXX, infoInvXX, Ncov;
	int D;
	double	oneD=1.0, moneD = -1.0,zeroD=0.0;
	double xybeta, *xvx;

	/// copy of LyLx
	Ncov = N*Nrep;
	F77_CALL(dcopy)(&Ncov,
			YXVinvYX, &oneI,
			copyLx, &oneI);

	// xvinvx submatrix
	xvx = &YXVinvYX[N*Nrep+Nrep];
	Ncov = N-Nrep;

	// invert X Vinv X
	// first cholesky

	//  cholesky X Vinv X
	F77_CALL(dpotrf)("L",
		&Ncov, xvx,
		&N, // Ncov by Ncov submatrix of N by N matrix
		&infoCholXX);

	*detXVinvX  = 0.0;
	for(D=0;D<Ncov;++D){
		*detXVinvX  += log(xvx[D*N+D]);
	}
	*detXVinvX *= 2;
// then invert
F77_NAME(dpotri)("L",
		&Ncov,
		xvx,
		&N,
		&infoInvXX);

// put beta hat in first rows (first Nrep cols still have LyLx)
// C= A B, A=xvx, B=LxLy
F77_NAME(dsymm)(
		"L", "L", // A on left, A in lower
		&Ncov, &Nrep, // C has Nrep rows, Ncov columns
		&oneD,
		xvx, &N,// XVinvX^(-1)
		&copyLx[Nrep],&N, // Lx
		&zeroD,
		&YXVinvYX[Nrep], &N //betahat, ldc
		);

//  blasBeta  C     alpha A     B
//     (1)  LyLy + (-1)  LxLy beta
F77_NAME(dgemm)(
		//	op(A), op(B),
		"T", "N",
		// nrows of op(A), ncol ob(B), ncol op(A) = nrow(opB)
  		&Nrep, &Nrep, &Ncov,
		// alpha
		&moneD,
		// A, lda(A), betahat
		&copyLx[Nrep], &N,
		// B, lda(B), LxLy
		&YXVinvYX[Nrep], &N, // betaHat
		// blasBeta
  		&oneD,
		// C, nrow(c)
		YXVinvYX, &N);
}


/*
 * logL given xisqTausq
// needs global variables
// Q, L, c, detTwo
// obsCovRot, Lx, YwkL, EwkL, DLx, YwkD, EwkD
// YXVYXglobal, YXYX, Nxy, Nobs,
 */

double logLoneNugget(double xisqTausq){

	double minusXisqTausq, zeroD=0.0, oneD=1.0;
	double *DYXVYX, result;
	int oneI=1, D;

	DYXVYX=YXVYXglobal;

	M_cholmod_factorize_p(
		Q,
		&xisqTausq, // beta
		(int*)NULL, 0 /*fsize*/,
		L, &c
	);


//Lx =
M_cholmod_solve2(
		CHOLMOD_L,
		L,
		obsCovRot,
		&Lx,
		&YwkL, &EwkL,
		&c);


// cross product
minusXisqTausq = -xisqTausq;

// - LxLx
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
		Lx->x, &Nobs,
		// beta
  		&zeroD,
		// C, nrow(c)
		DYXVYX, &Nxy);

// add the cross prod of data
F77_NAME(daxpy)(
		&Nxysq,
		&oneD,
		YXYX,
		&oneI,
		DYXVYX,
		&oneI);

detTwo[0] = M_chm_factor_ldetL2(L);

ssqFromXprod(
		DYXVYX, // Nxy by Nxy
		&detTwo[1],
		Nxy, Nrep, copyLx);

for(D=0;D<Nrep;++D){
	// using result as temporary variable
	result = log(DYXVYX[D*Nxy+D]);
// ml
logLtwo[D] = Nobs*result - Nobs*log(Nobs) + detTwo[0] - YrepAdd[D];
// reml
logLtwo[Nrep+D] = (Nobs-Ncov)*result -
		(Nobs-Ncov)*log(Nobs-Ncov) +
		detTwo[0] - detTwo[1] - YrepAdd[D];
}


//  now using DYXVYX as temporary variable
//if(Ltype){
//	DYXVYX = &logLtwo[Nxy];
//} else {
//	DYXVYX = &logLtwo[0];
//}

// find minimum element
//R_max_col(DYXVYX,&oneI, &Nrep,&D,&oneI);
//result = DYXVYX[D];

return result;

}

/*
 * callable function from R
 */
SEXP gmrfLik(
		SEXP QR,
		SEXP obsCovR,
		SEXP xisqTausq,
		SEXP YrepAddR
		){

	int DxisqTausq, NxisqTausq, Drep; // length(xisqTausq)
	double	oneD=1.0, zeroD=0.0;
	double *YXVYX, *determinant, *determinantForReml;
	double *m2logL, *m2logReL, *varHatMl, *varHatReml;
	SEXP resultR;
	CHM_DN obsCov;

	Ltype=0; // set to 1 for reml

	Nrep =LENGTH(YrepAddR);
	YrepAdd = REAL(YrepAddR);

	Nobs = INTEGER(GET_DIM(obsCovR))[0];
	Nxy = INTEGER(GET_DIM(obsCovR))[1];
	Ncov = Nxy - Nrep;
	NxisqTausq = LENGTH(xisqTausq);
	Nxysq = Nxy*Nxy;

	resultR = PROTECT(allocVector(REALSXP, Nxysq*NxisqTausq + 7*Nrep*NxisqTausq));


	YXVYX = REAL(resultR);
	determinant = &REAL(resultR)[Nxysq*NxisqTausq];
	determinantForReml = &REAL(resultR)[Nxysq*NxisqTausq + Nrep*NxisqTausq];
	m2logL = &REAL(resultR)[Nxysq*NxisqTausq + 2*Nrep*NxisqTausq];
	m2logReL = &REAL(resultR)[Nxysq*NxisqTausq + 3*Nrep*NxisqTausq];
	varHatMl = &REAL(resultR)[Nxysq*NxisqTausq + 4*Nrep*NxisqTausq];
	varHatReml = &REAL(resultR)[Nxysq*NxisqTausq + 5*Nrep*NxisqTausq];

	YXYX = (double *) calloc(Nxysq,sizeof(double));
	logLtwo = (double *) calloc(2*Nrep,sizeof(double));
	copyLx = (double *) calloc(Nxy*Nrep,sizeof(double));

	Q = AS_CHM_SP(QR);
	obsCov = AS_CHM_DN(obsCovR);
	M_R_cholmod_start(&c);

	// get some stuff ready

	// allocate Lx
	Lx = M_cholmod_copy_dense(obsCov,&c);

	// likelihood without nugget

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
			obsCov->x, &Nobs,
			// B, nrow(B)
	  		Lx->x, &Nobs,
			// beta
	  		&zeroD,
			// C, nrow(c)
			YXVYX, &Nxy);


	// Q = P' L D L' P
	L = M_cholmod_analyze(Q, &c);
	M_cholmod_factorize(Q,L, &c);

	// determinant
	determinant[0] = M_chm_factor_ldetL2(L);

	ssqFromXprod(
			YXVYX, // N by N
			determinantForReml,
			Nxy, Nrep,
			copyLx);

	for(Drep=0;Drep<Nrep;++Drep){
		determinant[Drep] = determinant[0];
		determinantForReml[Drep] = determinantForReml[0];

	m2logL[Drep] = Nobs*log(YXVYX[Drep*Nxy+Drep]) - Nobs*log(Nobs) -
			determinant[0] - YrepAdd[Drep];

	m2logReL[Drep] = (Nobs-Ncov)*log(YXVYX[Drep*Nxy+Drep]/(Nobs-Ncov)) +
			determinantForReml[0] - determinant[0] - YrepAdd[Drep];

	varHatMl[Drep] = YXVYX[Drep*Nxy+Drep]/Nobs;
	varHatReml[Drep] = YXVYX[Drep*Nxy+Drep]/(Nobs-Ncov);
	}
	// now with xisqTausq
	obsCovRot = M_cholmod_solve(CHOLMOD_P, L,obsCov,&c);

	// YXYX cross product of data
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



	for(DxisqTausq=1;DxisqTausq < NxisqTausq;++DxisqTausq){

		YXVYXglobal = &YXVYX[DxisqTausq*Nxysq];

		logLoneNugget(REAL(xisqTausq)[DxisqTausq]);

		// assign global values into their correct spot

		for(Drep=0;Drep<Nrep;++Drep){
			determinant[DxisqTausq*Nrep+Drep]=detTwo[0];
			determinantForReml[DxisqTausq*Nrep+Drep]=detTwo[1];

		m2logL[DxisqTausq*Nrep+Drep]  =  logLtwo[Drep] - determinant[0];
		m2logReL[DxisqTausq*Nrep+Drep] = logLtwo[Nrep+Drep] - determinant[0];

		varHatMl[DxisqTausq*Nrep + Drep] = YXVYXglobal[Drep*Nxy+Drep]/Nobs;
		varHatReml[DxisqTausq*Nrep + Drep] = YXVYXglobal[Drep*Nxy+Drep]/(Nobs-Ncov);
		}

	}

	M_cholmod_free_factor(&L, &c);
	M_cholmod_free_dense(&obsCovRot, &c);

	M_cholmod_free_dense(&Lx, &c);


// don't free Q because it's from an R object
//	M_cholmod_free_sparse(&Q, &c);

// don't free obsCov because it's from an R object
//	M_cholmod_free_dense(&obsCov, &c);

	free(copyLx);
	free(YXYX);
	free(logLtwo);
	M_cholmod_free_dense(&YwkL, &c);
	M_cholmod_free_dense(&YwkD, &c);
	M_cholmod_free_dense(&EwkL, &c);
	M_cholmod_free_dense(&EwkD, &c);

	M_cholmod_finish(&c);

	UNPROTECT(1);
	return resultR;
}


