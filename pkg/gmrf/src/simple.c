#include<R.h>
#include<R_ext/Applic.h>
#include<R_ext/Print.h>
#include<R_ext/Utils.h>


#include"Matrix.h"
#include"Matrix_stubs.c"
/** cholmod_common struct local to the lme4 package */
cholmod_common c;




SEXP simple(
		SEXP QR,
		SEXP obsCovR,
		SEXP xisqTausq
		){


	SEXP resultR;
	CHM_DN obsCov;
	CHM_SP Lmat;
	CHM_SP Q;
	CHM_FR L;


	resultR = PROTECT(allocVector(REALSXP, 4+20));


	Q = AS_CHM_SP(QR);
	obsCov = AS_CHM_DN(obsCovR);


	Rprintf("a");

	M_R_cholmod_start(&c);

	Rprintf("b");


	M_cholmod_finish(&c);
	Rprintf("c");



	UNPROTECT(1);
	return resultR;
}
