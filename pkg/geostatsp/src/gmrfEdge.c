#include"geostatsp.h"

SEXP gmrfEdge(
		SEXP QAA, // sparse symmetric
		SEXP QAB, // sparse
		SEXP outerDist, // dense symmetric
		SEXP params,
		SEXP gmrfParams
){

	SEXP result, typePrecision; // dense symmetric

	PROTECT(typePrecision = NEW_CHARACTER(1));
	SET_STRING_ELT(typePrecision, 0, mkChar("precision"));

	PROTECT(result = maternDistance(
    	outerDist,
		params,
    	typePrecision));

	// if length(gmrfParams) > 0
	// copy gmrfParams into QAA
	// copy gmrfParams into QAB
	// don't permute

	// cholInnerPrec =Cholesky(QAA,LDL=FALSE,perm=FALSE)
	// cholInnerPrec = Pt L Lt P
	// cholInnerPrec^(-1) = Pt L^(-1)t L^(-1) P
	// Aic = L^(-1) P QAB

	// if permute, reorder QAB
	/*
	Aic = solve(cholInnerPrec,
			QAB,
			system='P')
	 */

	/*
	Aic = solve(cholInnerPrec,
			QAB,
			system='L')
	 */


//	result = crossprod(Aic) + result
	// blas SYRK

	UNPROTECT(2);
	return result;
}
