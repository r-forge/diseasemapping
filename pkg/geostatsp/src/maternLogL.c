#include<R.h>
#include<Rmath.h>
#include<R_ext/Lapack.h>

void maternAniso(double *x, double *y, int *N,
		double *result,
		double  *range, double*shape, double *variance,
		double *anisoRatio, double *anisoAngleRadians,
		double *nugget, int *type
		);
void matern(double *distance, int *N,
		double *result,
		double *range, double *shape,
		double *variance,
		double *nugget, int *type);

void computeBoxCox(double *obsCov,
		// number of observations, number of datasets
		int *Nobs, int *Nrep,
		double *boxcox, // of length Nrep
		int *boxcoxType,
		// 0= do nothing
		// 1= do box-cox
		// 2= put log(y) in 2nd column of obscov
		//    ... and leave the first 2 columns as-is
		// 3= same as 2 but 2nd column already
		//    computed, as is sumLogL
		double *sumLogY
) {
//
	int D, Dbc, N, Nend;
	double *pLogY, *pRep, bcHere;

	if(!*boxcoxType){
		return;
	}

	N = *Nobs;

	if(*boxcoxType ==1){
		pLogY = obsCov; // logs go in first column
		Nend=-1;
	} else {
		pLogY = &obsCov[Nobs]; // logs go in 2nd col
		Nend = 1;
	}
	if(boxcoxType < 3){
		*sumLogY = 0.0;
		for(D=0;D<N;++D) {
			pLogY[D] = log(obscov[D]);
			*sumLogY += pLogY[D];
		}
	}

	Dbc = *Nrep-1;
	while(Dbc>Nend){
		pRep = &obsCov[Dbc*N];
		bcHere = boxcox[Dbc];
		for(D=0;D<N;++D) {
			pRep[D] = (exp(bcHere*pLogY[D]) - 1) /
					bcHere;
		}
	}

} // end box-cox


void maternLogLGivenChol(
		double *obsCov,
		int *Nobs, int *Nrep, int *Ncov,
		double *cholVariance,
		double *varMultHat, // a 2 by Nrep matrix (ml, reml
		double *betaHat, // an Ncov by Nrep matrix
		double *varBetaHat, // an Ncov by Ncov by Nrep array
		double *boxcox, // of length Nrep
		int *boxcoxType
		double *sumLogY, // of length Nrep
		double *logL, // a 4 by Nrep matrix
		  // ( ml fixed var, ml est v, reml f, reml est)
		) {

	int D, Ncol, infoCholLx, infoInvLx, NobsUse;
	double one, *pCov, *pObs;
	double *Lx, detLx, LxLy;

	one=1.0;

	if(*boxcoxType){
		computeBoxCox(obsCov,
				Nobs, Nrep,
				boxcox,
				boxcoxType,
				sumLogY);
	}
	if(*boxcoxType > 1)
		Ncol = *Nrep + *Ncov-2;
		pObs = &obsCov[*Nobs*2]];
		NrepUse = *Nrep - 2;
	} else {
		Ncol = *Nrep + *Ncov;
		pObs = obsCov;
		NrepUse = *Nrep;
	}

	Lx = (double *) calloc(*Ncov*(*Nobs),sizeof(double));
	LxLy = (double *) calloc(*Ncov*(*NrepUse),sizeof(double));

	// solve L x = b
	//      left, lower, not transposed, not diagonal
	F77_NAME(dtrsm)("L", "L", "N","N",
			Nobs, &Ncol,
			&one;
			cholVariance,
			Nobs,
	  		pObs,
			Nobs);

	//  cholCovInvXcross = Matrix::crossprod(cholCovInvX)
// transpose A, don't transpose B,
	pCov = &obsCov[*Nrep*(*Nobs)];
	F77_NAME(dgemm)('T', 'N',
			// nrows of op(A), ncol ob(B), ncol op(A) = nrow(opB)
	  		Ncov, Ncov, Nobs,
			&one,
	  		pCov, Nobs,
	  		pCov, Nobs,
	  		&one,
			Lx, Ncov);



	//cholCovInvXcrossInv =       Matrix::solve(cholCovInvXcross)
	//detCholCovInvXcross = Matrix::determinant(cholCovInvXcross)$modulus

	F77_CALL(dpotrf)("L", Ncov, Lx, Ncov, &infoCholLx);
	detLx=0;  // the log determinant
	for(D = 0; D < *Ncov; D++)
		detLx += log(Lx[D*(*Ncov)+D]);

	F77_NAME(dpotri)("L", Ncov,
			Lx, Ncov,
			&infoInvLx);

	//betaHat = as.vector(
	//      cholCovInvXcrossInv %*%
	//	 Matrix::crossprod(cholCovInvX, cholCovInvY))
	F77_NAME(dgemm)('T', 'N',
			// nrows of op(A), ncol ob(B), ncol op(A) = nrow(opB)
	  		Ncov, NrepUse, Nobs,
			&one,
	  		pCov, Ncov,
	  		pObs, Nobs,
	  		&one,
			LxLy, Ncov);
	F77_NAME(dgemm)('T', 'N',
			// nrows of op(A), ncol ob(B), ncol op(A) = nrow(opB)
	  		Ncov, Ncov, Nobs,
			&one,
	  		Lx, Nobs,
			LxLy, Nobs,
	  		&one,
			betaHat, Ncov);



	  free(Lx);
	  free(LxLy);

}

void maternLogL(
		double *obsCov, int *Nobs, int *Ny, int *Ncov,
		double *xcoord, double *ycoord,
		// either xcoord is distances, length Nobs*Nobs
		// and ycoord is ignored, or xcoord and ycoord
		// are vecctors of length Nobs
		double *range, double *shape,
		double *variance,
		double *anisoRatio,
		double *anisoAngleRadians,
		double *nugget,
		double *boxcox, // of length Ny
		double *sumLogY, // of length Ny
		double *logL, // a 4 by Ny matrix
		  // ( ml fixed var, ml est v, reml f, reml est)
		double *varMultHat, // a 2 by Ny matrix (ml, reml
		double *betaHat, // an Ncov by Ny matrix
		double *varBetaHat, // an Ncov by Ncov by Ny array
		int *useBoxCox,
		int *aniso
		) {

  char *lower = "L";
  int n = *Nobs , i;
  int nn = n*n;
  int info, two;
  double *corMat, logDet;

  two=2;
  corMat = (double *) calloc(nn,sizeof(double));

  if(*aniso) {
	  maternAniso(xcoord,ycoord,Nobs,corMat,range,shape,
			  variance,anisoRatio,
			  anisoAngleRadians,nugget,&two);
  } else {
	  matern(xcoord,Nobs,corMat,
			range,shape,
			variance,nugget,&two);
  }
  logDet = *shape;

  // solve L x = b
  /*
  F77_NAME(dtrsm)(const char *side, const char *uplo,
  		const char *transa, const char *diag,
  		const int *m, const int *n, const double *alpha,
  		const double *a, const int *lda,
  		double *b, const int *ldb);

  F77_NAME(dgemm)(const char *transa, const char *transb, const int *m,
  		const int *n, const int *k, const double *alpha,
  		const double *a, const int *lda,
  		const double *b, const int *ldb,
  		const double *beta, double *c, const int *ldc);

  F77_NAME(dpotrf)(const char* uplo, const int* n,
  		 double* a, const int* lda, int* info);
  		 */
  /* DPOTRI - compute the inverse of a real symmetric positive */
  /* definite matrix A using the Cholesky factorization A = U**T*U */
  /* or A = L*L**T computed by DPOTRF */

//  F77_NAME(dpotri)(const char* uplo, const int* n,
//  		 double* a, const int* lda, int* info);


}
