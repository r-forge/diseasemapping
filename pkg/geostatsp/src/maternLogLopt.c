#include"geostatsp.h"
#include<R_ext/Applic.h>
#include<R_ext/Print.h>

int *SparamOpt, *limTypeOpt; // number of params, followed by indices of params to be optimized
double *paramOpt;
// 0 nugget, 1 variance,
// 2 range, 3 shape,
// 4 anisoRatio, 5 ansioAngleRadians,
// 6 boxcox
double *lower, *upper, *parscale, *ndeps;

double *obsCovOpt, *obsCovCopy; // pointer to the data
double *obsForBoxcoxOpt; // pointer three column matrix of Y, logY and boxcox Y

const double *xcoordOpt, *ycoordOpt;
double *corMatOpt;

double boxcoxParamOpt[9] = {1,0,1,-9,-9,-9,-9,-9,-9};
int anisoOpt, LtypeOpt, boxcoxTypeOpt;
int Nopt[3];// Nobs, 1, Ncov
int NboxcoxOpt[2], doBoxcoxOpt; // Nobs, 3

double *totalVarHatOpt;
int NforBoxCoxOpt[3];
double  *betaHatOpt, *varBetaHatOpt;
double *LxLyOpt;


double maternLogLObj(
		int junk,
		double *paramArg,
		void *ex
		) {
	int zero=0, oneI=1, three=3, maternType = 2, Dparam;
	int NforCopy;
	double determinants[2], logL[4];
	// copy param to fullParam

	NforCopy = Nopt[0]*(Nopt[1]+Nopt[2]);

	for(Dparam=0;Dparam<SparamOpt[0];++Dparam){
		paramOpt[SparamOpt[Dparam+1]]=paramArg[Dparam];
	}

	if(doBoxcoxOpt){
		boxcoxParamOpt[2] = paramOpt[6];
		computeBoxCox(
			obsForBoxcoxOpt,
			NboxcoxOpt,
			boxcoxParamOpt,
			&three
			);
	}

	// make a copy of the data
	F77_NAME(dcopy)(&NforCopy,
			obsCovOpt, &oneI,
			obsCovCopy, &oneI);

	maternForL(
		xcoordOpt, ycoordOpt,
		Nopt,corMatOpt,
		paramOpt,
		&anisoOpt,
		&zero,// don't ignore nugget
		&maternType,//chol of variance matrix
		determinants);

	maternLogLGivenChol(
		obsCovCopy,
		Nopt,  // Nobs, 1, Ncov,
		corMatOpt,
		logL, // a 1 by Nrep matrix
		betaHatOpt, // an Ncov by Nrep matrix
		varBetaHatOpt, // an Ncov by Ncov by Nrep array
		determinants, // detVarHalf, detCholCovInvXcrossHalf
		LxLyOpt);

	logLfromComponents(
				Nopt,
				boxcoxParamOpt,
				&doBoxcoxOpt,
				logL,
				totalVarHatOpt,
				determinants,
				&LtypeOpt
		);
//	Rprintf("\n p ");
//	for(Dparam=0;Dparam<SparamOpt[0];++Dparam)
//		Rprintf(" %f ", paramArg[Dparam] );
//	Rprintf("\n pf ");
//	for(Dparam=0;Dparam<7;++Dparam)
//		Rprintf(" %f ", paramOpt[Dparam] );
//	Rprintf("l %f ", logL[0]);

	return logL[0];
}

void maternLogLgr(
		int junk,
		double *paramArg,
		double *result,
		void *ex
) {
	int Dpar, oneI=1, Nparam;
	double *parHere, deltaPar, *fullGr;

	Nparam = SparamOpt[0];

	parHere = (double *) calloc(Nparam,sizeof(double));
	fullGr = (double *) calloc(6*Nparam,sizeof(double));


	for(Dpar=0;Dpar < Nparam;++Dpar){
//		Rprintf(" %d %d ", Dpar, Nparam);
		deltaPar = parscale[Dpar]* ndeps[Dpar];

		//lower
		F77_NAME(dcopy)(&Nparam,
				paramArg, &oneI,
				parHere, &oneI);

		parHere[Dpar]= paramArg[Dpar]-deltaPar;
		if(limTypeOpt[Dpar] == 1 | limTypeOpt[Dpar] == 2){
			parHere[Dpar]=fmax(
					parHere[Dpar],
					lower[Dpar]
			);
		}
		fullGr[Dpar+1*Nparam]=parHere[Dpar];

		fullGr[Dpar+2*Nparam] = maternLogLObj(
				junk,parHere, ex);

		//upper
		parHere[Dpar]=paramArg[Dpar]+deltaPar;
		if(limTypeOpt[Dpar] == 3 | limTypeOpt[Dpar] == 2){
			parHere[Dpar] = fmin(
					parHere[Dpar],
					upper[Dpar]
			);
		}
		fullGr[Dpar+3*Nparam]=parHere[Dpar];

		fullGr[Dpar+4*Nparam] = maternLogLObj(
				junk,parHere, ex);
		fullGr[Dpar+5*Nparam] =fullGr[Dpar+4*Nparam] -
				fullGr[Dpar+2*Nparam];
		result[Dpar] = fullGr[Dpar+5*Nparam]/
				(fullGr[Dpar+3*Nparam] - fullGr[Dpar+1*Nparam]);

		/*
  		Rprintf(" l%f u%f d%f ",
				fullGr[Dpar+2*Nparam],
				fullGr[Dpar+4*Nparam],
				result[Dpar] );
		*/
	}

//	Rprintf("  \n");

	free(parHere);
	free(fullGr);

}


void maternLogLOpt(
		double *fullParam,// nugget, variance,
        // range, shape,
        // anisoRatio, ansioAngleRadians, boxcox
		int *Sparam,
		// vector of 0 and 1, depending on whether corresponding parameter
		// in fullParam is to be optimized
		double *obsCov,
		const double *xcoord,
		const double *ycoord,
		const int *aniso,
		const int *N,// Nobs, Nrep, Ncov
		int *Ltype,
		int *scalarsInt,
		double *scalarsF,
		double *parLim,
		int *limType,
		char **msg
		// 0=ml, var estimated
		// 1=reml, var estimated
		//2=ml, var fixed
		// 3=reml, var fixed
		// on exit, info from chol of matern
) {

	double *paramArg, result;
	int two=2, oneI=1,Dparam, DparamForOpt, optimFail;
	int fncount, grcount, junk;
	void *optimEx;
	char themsg[100];

	// assign the global variables (which end in Opt)
	xcoordOpt=xcoord;
	ycoordOpt=ycoord;

	Nopt[0]	= N[0];
	Nopt[1]	= 1;
	Nopt[2]	= N[2];

	paramOpt=fullParam;
	anisoOpt=*aniso;


	SparamOpt = (int *) calloc(7,sizeof(int));

	// create a vector of indices of parameters to estimate
	DparamForOpt = 1;
	for(Dparam=0;Dparam<6;++Dparam){
		if(Sparam[Dparam]){
			SparamOpt[DparamForOpt] = Dparam;
			++DparamForOpt;
		}
	}
	// first entry is number of parameters to be estimated
	SparamOpt[0]=DparamForOpt-1;

	paramArg = (double *) calloc(SparamOpt[0],sizeof(double));
	totalVarHatOpt = (double *) calloc(4,sizeof(double));
	betaHatOpt= (double *) calloc(N[2],sizeof(double));
	varBetaHatOpt = (double *) calloc(N[2]*N[2],sizeof(double));
	LtypeOpt = *Ltype;


	// create vector of values for parameters to be estimated
	for(Dparam=0;Dparam<SparamOpt[0];++Dparam){
		paramArg[Dparam] = fullParam[SparamOpt[Dparam+1]];
	}

	// for optim
	lower = parLim;
	upper = &parLim[SparamOpt[0]];
	parscale = &parLim[SparamOpt[0]*2];
	ndeps = &parLim[SparamOpt[0]*3];
	limTypeOpt = limType;


	// prepare boxcox stuff
	// first two columns of obscov are Y and log(Y)
	// likelihood not computed for them.
	obsCovOpt = &obsCov[ 2*N[0] ];
	// but they are used for computing box-cox transform
	obsForBoxcoxOpt = obsCov;
	// are we doing box-cox?
	if(Sparam[6]){ // yes, it's being optimized
		doBoxcoxOpt=1;
		// put log(y) in position 2
		NboxcoxOpt[0] = N[0];
		NboxcoxOpt[1] = 3;
		boxcoxTypeOpt=2;
		computeBoxCox(
				obsCov,
				NboxcoxOpt,
				boxcoxParamOpt,
				&boxcoxTypeOpt
		);
		boxcoxTypeOpt=3;
	} else { // not being optimized
		doBoxcoxOpt=0;
		// but check if it's not 1
		if(fabs(fullParam[6]-1)>0.001){
			NboxcoxOpt[0] = N[0];
			NboxcoxOpt[1] = 3;
			boxcoxTypeOpt=1;
			computeBoxCox(
				obsCov,
				NboxcoxOpt,
				boxcoxParamOpt,
				&boxcoxTypeOpt
				);
			boxcoxTypeOpt=4;
		} else {
			boxcoxTypeOpt=0;
			// put observations in column 3
			F77_NAME(dcopy)(N,
					obsCov, &oneI,
					obsCovOpt, &oneI);
			boxcoxParamOpt[8]=0.0;
		}
	}


// allocate memory

	obsCovCopy = (double *) calloc(N[0]*(1+N[2]),sizeof(double));
	// enough memory for covariates and one vector of observations
	corMatOpt = (double *) calloc(N[0]*N[0],sizeof(double));
	LxLyOpt = (double *) calloc(100+N[1]*N[2],sizeof(double));

#ifdef UNDEF
	resultGr = (double *) calloc(SparamOpt[0]*6,sizeof(double));

	result = maternLogLObj(junk,paramArg, optimEx);
	totalVarHat[3] = result;
	totalVarHat[1] =maternLogLObj(junk,paramArg, optimEx);

	maternLogLgr(
			junk,
			paramArg,
			resultGr,
			optimEx);

	for(Dparam=0;Dparam<SparamOpt[0]*6;++Dparam){
		parLim[Dparam] = resultGr[Dparam];
	}
# endif

	lbfgsb(
			SparamOpt[0],//int n,
			scalarsInt[4],//int lmm,
			paramArg,//double *x,
			lower,
		    upper,
			limType,//int *nbd,
			&result,//double *Fmin,
			maternLogLObj,//optimfn fn,
			maternLogLgr,
			&optimFail,//int *fail,
			optimEx,//void *ex,
			scalarsF[6],//double factr,
			scalarsF[7],//double pgtol,
			&fncount,
			&grcount,
			scalarsInt[1],//int maxit,
			themsg,//char *msg,
			scalarsInt[0],//int trace,
			scalarsInt[2]//int nREP
		  );

	strcpy(*msg, themsg);

	scalarsInt[0] = optimFail;
	scalarsInt[1] = fncount;
	scalarsInt[2] = grcount;
	scalarsF[0] = result;
	scalarsF[1] = totalVarHatOpt[0];


	for(Dparam=0;Dparam<SparamOpt[0];++Dparam){
		paramOpt[SparamOpt[Dparam+1]]=paramArg[Dparam];
	}

	// make a copy of betahat
	F77_NAME(dcopy)(&N[2],
			betaHatOpt, &oneI,
			parLim, &oneI);
	// and varBetaHat
	Dparam = N[2]*N[2];
	F77_NAME(dcopy)(&Dparam,
			varBetaHatOpt, &oneI,
			&parLim[N[2]], &oneI);

	free(corMatOpt);
	free(obsCovCopy);
	free(LxLyOpt);
	free(paramArg);
	free(SparamOpt);
	free(totalVarHatOpt);
	free(betaHatOpt);
	free(varBetaHatOpt);
}

