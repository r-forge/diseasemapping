/*
 * needs modified matrix package
 * svn checkout svn://svn.r-forge.r-project.org/svnroot/matrix/pkg/Matrix Matrix
 R CMD build --no-build-vignettes Matrix
 no vignettes  because I don't have texi2pdf
 */

#include"geostatsp.h"
//#include<R.h>
//#include<Rmath.h>
//#include<R_ext/Lapack.h>
//#include<R_ext/Applic.h>
//#include<R_ext/Print.h>
//#include<R_ext/Utils.h>
//#include<R_ext/Rdynload.h>
//#include<Matrix.h> it's in matrix_stubs.c
//#include<Matrix_stubs.c>

double Brent_fmin(double ax, double bx, double (*f)(double, void *),
		void *info, double tol);


/*
	a function for interfacing to the Matrix package
	hopefully this will be in Matrix_stubs.c soon
 */


/*
 *  global variables, so that an optimizer can be build eventually
 */
CHM_SP Q;
CHM_FR L;
CHM_DN obsCovRot, Lx;
CHM_DN YwkL, EwkL, YwkD, EwkD; // workspaces
cholmod_common c;
double *YXVYXglobal,  *YXYX, *YrepAdd;
double *copyLx;
int Nxy, Nobs, Ncov, Nrep, Nxysq;
int Ltype;
int DxisqTausq, NxisqTausq;

/*
 * compute sums of squares from cross products
 *
 */
void ssqFromXprod(
		double *YXVinvYX, // N by N
		double *detXVinvX,
		const int N, // Nxy
		const int Nrep,
		double *copyLx
){

	int oneI=1, infoCholXX, infoInvXX, Ncov;
	int D;
	double	oneD=1.0, moneD = -1.0,zeroD=0.0;
	double *xvx;

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

	// put beta hat in first column (first Nrep rows still have LyLx)
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




double logLoneNuggetMoreArguments(
		double xisqTausq,
		double *DYXVYX,
		double *determinant,
		double *determinantForReml,
		double *m2logL,
		double *m2logReL
){
	// Q, L, Lx, copyLx, c, obsCovRot, YXYX, YrepAdd are global
	// and int Nxy, Nobs, Nxysq, Nrep;
	// as are YwkL, EwkL, YwkD, EwkD; // workspaces

	double minusXisqTausq, zeroD=0.0, oneD=1.0;
	double result;
	int oneI=1, D;

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

	*determinant = M_chm_factor_ldetL2(L);

	ssqFromXprod(
			DYXVYX, // Nxy by Nxy
			determinantForReml,
			Nxy, Nrep, copyLx);

	for(D=0;D<Nrep;++D){
		// using result as temporary variable
		result = log(DYXVYX[D*Nxy+D]);
		// ml
		m2logL[D] = Nobs*result - Nobs*log(Nobs) +
				*determinant - YrepAdd[D];
		// reml
		m2logReL[D] = (Nobs-Ncov)*result -
				(Nobs-Ncov)*log(Nobs-Ncov) +
				*determinant - *determinantForReml - YrepAdd[D];
	}


	//  now using DYXVYX as temporary variable
	if(Ltype){
		DYXVYX = m2logReL;
	} else {
		DYXVYX = m2logL;
	}

	// find minimum element
	R_max_col(DYXVYX,&oneI, &Nrep,&D,&oneI);
	result = DYXVYX[D];

	return result;

}

/*
 * logL given xisqTausq
// needs global variables
// Q, L, c, detTwo
// obsCovRot, Lx, YwkL, EwkL, DLx, YwkD, EwkD
// YXVYXglobal, YXYX, Nxy, Nobs,
 */

double logLoneNugget(double xisqTausq, void *nothing){

	return logLoneNuggetMoreArguments(
			xisqTausq,
			&YXVYXglobal[DxisqTausq*Nxysq], // data cross prod
			&YXVYXglobal[Nxysq*NxisqTausq+
						 DxisqTausq*Nrep], // det
						 &YXVYXglobal[Nxysq*NxisqTausq +
									  Nrep*NxisqTausq+
									  DxisqTausq*Nrep], // reml det
									  &YXVYXglobal[Nxysq*NxisqTausq +
												   2*Nrep*NxisqTausq+
												   DxisqTausq*Nrep], // logL
												   &YXVYXglobal[Nxysq*NxisqTausq +
																3*Nrep*NxisqTausq+
																DxisqTausq*Nrep] // reLogL
	);

}

// logL for calling form Brent_fmin
// this needs YXVYXstart as well as the previous
// function's global arguments
double logLoneLogNugget(double logXisqTausq, void* nothing){

	double result;

	if(DxisqTausq>=NxisqTausq){
		DxisqTausq=1;
	}

	result = logLoneNugget(exp(logXisqTausq), nothing);
	DxisqTausq++;

	return(result);
}

/*
 * callable function from R
 */
SEXP gmrfLik(
		SEXP QR,
		SEXP obsCovR,
		SEXP xisqTausq, // first element must be zero
		SEXP reml,
		SEXP YrepAddR,// Jacobian, length used for Nrep
		SEXP optParam  // set first element negative for verbose
){

	int Drep, dooptim, verbose=0; // length(xisqTausq)
	double	oneD=1.0, zeroD=0.0;
	double optTol, optMin, optMax; // default interval for optimizer
	int NxisqMax = 100; // number of xisqTausq's to retain when optimizing
	double *YXVYX, *determinant, *determinantForReml;
	double *m2logL, *m2logReL, *varHatMl, *varHatReml, *resultXisqTausq;
	double nothing=0.0;
	SEXP resultR;
	CHM_DN obsCov;

	Ltype=*INTEGER(reml); // set to 1 for reml

	// default interval for optimizer
	optMin = REAL(optParam)[0];
	optMax = REAL(optParam)[1];
	optTol = REAL(optParam)[2];// tolerence for fmin_bren

	Nrep =LENGTH(YrepAddR);
	YrepAdd = REAL(YrepAddR);

	Nobs = INTEGER(GET_DIM(obsCovR))[0];
	Nxy = INTEGER(GET_DIM(obsCovR))[1];
	Ncov = Nxy - Nrep;
	Nxysq = Nxy*Nxy;

	NxisqTausq = LENGTH(xisqTausq);
	// if length zero, do optimization
	dooptim=!NxisqTausq;

	if(dooptim){
		NxisqTausq = NxisqMax;
	} else {
		if(REAL(xisqTausq)[0]>0.0) {
			warning(
					"first element of xisqTausq is %f, must be zero\n",
					REAL(xisqTausq));
		}
		if(optMin<0.0) {
			verbose = 1;
			Rprintf("d Nxixqtauxq %d, Nobs %d, Nrep %d, Ncov %d\n",
					NxisqTausq, Nobs, Nrep, Ncov);
		}
	}

	resultR = PROTECT(allocVector(REALSXP,
			Nxysq*NxisqTausq + 8*Nrep*NxisqTausq));

	for(Drep = 0; Drep < LENGTH(resultR);++Drep){
		REAL(resultR)[Drep] = NA_REAL;
	}

	/*
	 * if NxisqTausq=2, Ny=1, Nx=1, Nxysq=4,
	 *  result has length 2*2*2+8*2=24
	 *  result[0-3] = t(obscov) Q obscov
	 *  result[4-7] = ?
	 *  result[8] = determinant
	 *  result[10] = detReml
	 *  result[12] = m2logL
	 *  result[14] = m2logLreml
	 *  result[16] = varhatMl
	 *  result[18] = varhatReml
	 *  result[20] = resultxisqTausq
	 */
	YXVYX = &REAL(resultR)[0];
	YXVYXglobal =  &REAL(resultR)[0];//YXVYX;

	determinant = &REAL(resultR)[Nxysq*NxisqTausq];
	determinantForReml = &REAL(resultR)[Nxysq*NxisqTausq + Nrep*NxisqTausq];
	m2logL = &REAL(resultR)[Nxysq*NxisqTausq + 2*Nrep*NxisqTausq];
	m2logReL = &REAL(resultR)[Nxysq*NxisqTausq + 3*Nrep*NxisqTausq];
	varHatMl = &REAL(resultR)[Nxysq*NxisqTausq + 4*Nrep*NxisqTausq];
	varHatReml = &REAL(resultR)[Nxysq*NxisqTausq + 5*Nrep*NxisqTausq];
	resultXisqTausq = &REAL(resultR)[Nxysq*NxisqTausq + 6*Nrep*NxisqTausq];

	YXYX = (double *) calloc(Nxysq,sizeof(double));
	copyLx = (double *) calloc(Nxy*Nrep,sizeof(double));


	Q = AS_CHM_SP(QR);
	obsCov = AS_CHM_DN(obsCovR);
	M_R_cholmod_start(&c);

	// get some stuff ready
	if(verbose) {
		Rprintf("starting\n");
	}
	/* .. now allocate Lx .. */
	Lx = M_cholmod_copy_dense(obsCov,&c);


	// for YX Vinv YX
	// Lx is Q YX
	M_cholmod_sdmult(
			Q,
			0, &oneD, &zeroD, // transpose, scale, scale
			obsCov,Lx,// in, out
			&c);
	// now compute t(YX) Lx
	if(verbose) {
		Rprintf("put t(obscov) Q obscov in result\n");
	}
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
	if(verbose) {
		Rprintf("first cholesky\n");
	}
	// Q = P' L D L' P
	L = M_cholmod_analyze(Q, &c);
	M_cholmod_factorize(Q,L, &c);

	// determinant
	determinant[0] = M_chm_factor_ldetL2(L);
	resultXisqTausq[0]= NA_REAL;

	if(verbose) {
		Rprintf("likelihood without nugget, det %f\n",
				determinant[0]);
	}

	ssqFromXprod(
			YXVYX, // N by N
			determinantForReml,
			Nxy, Nrep,
			copyLx);
#ifdef UNDEF
	for(Drep=0;Drep<Nrep;++Drep){
		determinant[Drep] = determinant[0];
		determinantForReml[Drep] = determinantForReml[0];

		m2logL[Drep] = Nobs*log(YXVYX[Drep*Nxy+Drep]) -
				Nobs*log(Nobs) -
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



	//	Rprintf("done zero %d\n", dooptim);


	if(dooptim){
		//		Rprintf(" opt %f %f %f\n", optMin, optMax, optTol);

		DxisqTausq=1;
		// do optimizer
		Brent_fmin(
				optMin, optMax,
				logLoneLogNugget,
				&nothing, optTol);



		NxisqTausq = DxisqTausq;

	} 	else { // not optimizing

		for(DxisqTausq=1;DxisqTausq < NxisqTausq;++DxisqTausq){
			if(verbose) {
				Rprintf("DxisqTausq %d, xisqTausq %f\n",
						DxisqTausq, REAL(xisqTausq)[DxisqTausq]);
			}
			logLoneNuggetMoreArguments(
					REAL(xisqTausq)[DxisqTausq],
					&YXVYX[DxisqTausq*Nxysq], // data cross prod
					&YXVYX[Nxysq*NxisqTausq+
						   DxisqTausq*Nrep], // det
					&YXVYX[Nxysq*NxisqTausq +
						  Nrep*NxisqTausq+
						  DxisqTausq*Nrep], // reml det
					&YXVYX[Nxysq*NxisqTausq +
						 2*Nrep*NxisqTausq+
						 DxisqTausq*Nrep], // logL
					&YXVYX[Nxysq*NxisqTausq +
						3*Nrep*NxisqTausq+
						DxisqTausq*Nrep] // reLogL
			);

		}
	}
	// assign global values into their correct spot

	for(DxisqTausq=0;DxisqTausq < NxisqTausq;++DxisqTausq){
		for(Drep=0;Drep<Nrep;++Drep){
			resultXisqTausq[DxisqTausq*Nrep + Drep] =
					REAL(xisqTausq)[DxisqTausq];
		}
	}
	for(DxisqTausq=1;DxisqTausq < NxisqTausq;++DxisqTausq){

		for(Drep=0;Drep<Nrep;++Drep){

			/*
			determinant[DxisqTausq*Nrep+Drep] =
					determinant[DxisqTausq*Nrep];

			determinantForReml[DxisqTausq*Nrep+Drep]=
					determinantForReml[DxisqTausq*Nrep];

			m2logL[DxisqTausq*Nrep+Drep]  -=
					determinant[0];

			m2logReL[DxisqTausq*Nrep+Drep] -=
					determinant[0];
		 */
			varHatMl[DxisqTausq*Nrep + Drep] =
					YXVYX[DxisqTausq*Nxysq+Drep*Nxy+Drep]/Nobs;

			varHatReml[DxisqTausq*Nrep + Drep] =
					YXVYX[DxisqTausq*Nxysq+Drep*Nxy+Drep]/(Nobs-Ncov);

			}

	}
#endif

	M_cholmod_free_factor(&L, &c);
	M_cholmod_free_dense(&obsCovRot, &c);

	M_cholmod_free_dense(&Lx, &c);


	// don't free Q because it's from an R object
	//	M_cholmod_free_sparse(&Q, &c);

	// don't free obsCov because it's from an R object
	//	M_cholmod_free_dense(&obsCov, &c);

	free(copyLx);
	free(YXYX);
	M_cholmod_free_dense(&YwkL, &c);
	M_cholmod_free_dense(&YwkD, &c);
	M_cholmod_free_dense(&EwkL, &c);
	M_cholmod_free_dense(&EwkD, &c);

	M_cholmod_finish(&c);

	UNPROTECT(1);
	return resultR;
}


