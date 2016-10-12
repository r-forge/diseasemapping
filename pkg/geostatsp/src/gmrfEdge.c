#include"geostatsp.h"

SEXP gmrfEdge(
		SEXP QAA, // sparse symmetric
		SEXP QAB, // sparse
		SEXP outerDist, // dense symmetric
		SEXP params
){

	SEXP result, typePrecision; // dense symmetric

	PROTECT(typePrecision = NEW_CHARACTER(1));
	SET_STRING_ELT(typePrecision, 0, mkChar("precision"));

	PROTECT(result = maternDistance(
    	outerDist,
		params,
    	typePrecision));

	//	cholInnerPrec =Cholesky(QAA,LDL=FALSE,perm=FALSE)

	/*
	Aic = solve(cholInnerPrec,
			QAB,
			system='L')
	 */

//	AQinvA = crossprod(Aic)
//	result = result + AQinvA

	UNPROTECT(2);
	return result;
}
