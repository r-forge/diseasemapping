#include"geostatsp.h"
#include<R_ext/Applic.h>

int SparamOpt[8], *limTypeOpt; // number of params, followed by indices of params to be optimized
double *paramOpt;
// 0 nugget, 1 variance,
// 2 range, 3 shape,
// 4 anisoRatio, 5 ansioAngleRadians,
// 6 boxcox
double *lower, *upper, *parscale, *ndeps;

double *obsCovOpt; // pointer to the data
double *obsForBoxcoxOpt; // pointer three column matrix of Y, logY and boxcox Y

const double *xcoordOpt, *ycoordOpt;
double *corMatOpt;

double boxcoxParamOpt[9] = {1,0,1,0,0,0};
int anisoOpt, LtypeOpt, boxcoxTypeOpt;
int Nopt[3];// Nobs, 1, Ncov
int NboxcoxOpt[2], doBoxcoxOpt; // Nobs, 3

double *totalVarHatOpt;
int NforBoxCoxOpt[3];
double  *betaHatOpt, *varBetaHatOpt;
double *LxLyOpt;


double maternLogLObj(
		int *junk,
		double *paramArg,
		void *ex
		) {
	int three=3, maternType = 2, zero=0, Dparam;
	double determinants[2], logL[4];
	// copy param to fullParam

	for(Dparam=0;Dparam<*SparamOpt;++Dparam){
		paramOpt[SparamOpt[Dparam]]=paramArg[Dparam];
	}

	if(doBoxcoxOpt){
		boxcoxParamOpt[2] = paramOpt[6];
		computeBoxCox(
				obsCovOpt,
			NboxcoxOpt,
			boxcoxParamOpt,
			&three
			);
	}

	maternForL(
		xcoordOpt, ycoordOpt,
		Nopt,corMatOpt,
		paramOpt,
		&anisoOpt,
		&zero,// don't ignore nugget
		&maternType,//chol of variance matrix
		determinants);

	maternLogLGivenChol(
		obsCovOpt,
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
	return logL[0];
}

double maternLogLgr(
		int *junk,
		double *paramArg,
		double *result,
		void *ex
) {
	int Dpar, oneI=1, Nparam, NparamFull, DparMat;
	double *parHere, deltaPar, deltaHere;

	Nparam = SparamOpt[0];
	NparamFull = Nparam; // have decided only to include pars to be optimized

	parHere = (double *) calloc(Nparam,sizeof(double));

// fill in parMat


	for(Dpar=0;Dpar < Nparam;++Dpar){
		deltaPar = parscale[Dpar]* ndeps[Dpar];

		//lower
		F77_NAME(dcopy)(Nparam,
				paramArg, &oneI,
				&parHere, &oneI);

		parHere[Dpar]= paramArg[Dpar]-deltaPar;
		if(limTypeOpt == 1 | limTypeOpt == 2){
			parHere[Dpar]=fmax(
					parHere[Dpar],
					lower[Dpar]
			);
		}

		deltaHere  = - parHere[Dpar];
		result[Dpar] = -maternLogLObj(
				junk,parHere, ex);

		//upper
		parHere[Dpar]=paramArg[Dpar]+deltaPar;
		if(limTypeOpt == 3 | limTypeOpt == 2){
			parMat[DparMat] = fmin(
					parMat[DparMat],
					upper[Dpar]
			);
		}

		deltaHere  += parHere[Dpar];
		result[Dpar] += maternLogLObj(
				junk,parHere, ex);

		result[Dpar] = result[Dpar]/deltaHere;
	}

	free(parHere);
 }


double maternLogLOpt(
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
		double *totalVarHat,
		double *betaHat,
		double *varBetaHat,
		int *Ltype,
		int *scalarsInt,
		double *scalarsF,
		double *parLim,
		int *limType,
		char *msg;
		// 0=ml, var estimated
		// 1=reml, var estimated
		//2=ml, var fixed
		// 3=reml, var fixed
		// on exit, info from chol of matern
) {

	double paramArg[7], result;
	int two=2, Dparam, DparamForOpt, optimFail, fncount, grcount, *junk;
	void *optimEx;

	// assign the global variables (which end in Opt)
	xcoordOpt=xcoord;
	ycoordOpt=ycoord;

	Nopt[0]	= N[0];
	Nopt[1]	= 1;
	Nopt[2]	=N[2];

	// for optim
	lower = parLim;
	upper = &parLim[Sparam[0]];
	parscale = &parLim[Sparam[0]*2];
	ndeps = &parLim[Sparam[0]*3];
	limTypeOpt = limType;

	paramOpt=fullParam;
	anisoOpt=*aniso;

	totalVarHatOpt = totalVarHat;
	betaHatOpt=betaHat;
	varBetaHatOpt = varBetaHat;
	LtypeOpt = *Ltype;

	// create a vector of indices of parameters to estimate
	DparamForOpt = 1;
	for(Dparam=0;Dparam<7;++Dparam){
		if(Sparam[Dparam]){
			SparamOpt[DparamForOpt] = Dparam;
			++DparamForOpt;
		}
	}
	// first entry is number of parameters to be estimated
	SparamOpt[0]=DparamForOpt-1;

	// create vector of values for parameters to be estimated
	for(Dparam=0;Dparam<*SparamOpt;++Dparam){
		paramArg[Dparam] = fullParam[SparamOpt[Dparam]];
	}


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
		NboxcoxOpt[1] = 2;
		boxcoxTypeOpt=2;
		computeBoxCox(
				obsForBoxcoxOpt,
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
			NboxcoxOpt[1] = 1;
			boxcoxTypeOpt=1;
			computeBoxCox(
				obsCovOpt,
				NboxcoxOpt,
				boxcoxParamOpt,
				&boxcoxTypeOpt
				);
			boxcoxTypeOpt=4;
		} else {
			boxcoxTypeOpt=0;
		}
	}


// allocate memory
	corMatOpt = (double *) calloc(N[0]*N[0],sizeof(double));
	LxLyOpt = (double *) calloc(N[0]*(N[1]),sizeof(double));

//	result = maternLogLObj(junk,paramArg, optimEx);
//	totalVarHat[0] = result;

	lbfgsb(
			Sparam[0],//int n,
			scalarInt[4],//int lmm,
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
			scalarInt[1],//int maxit,
			msg,//char *msg,
			scalarInt[0],//int trace,
			scalarInt[2]//int nREP
		  );

	scalarsInt[0] = *optimFail;
	scalarsInt[1] = *fncount;
	scalarsInt[2] = *grcount;

	free(corMatOpt);
	free(LxLyOpt);
}

