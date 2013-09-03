#include<R.h>
#include<Rmath.h>


void matern(double *distance, long *N,
		double *range, double *rough, double *variance) {

	int D, N2;
	double xscale, varscale,  thisx;

    long nb,  ize;
    double *bk, alpha;

	ize = 1L;
	alpha = *rough;

// code stolen from R's src/nmath/bessel_k.c
	nb = 1+ (long)floor(alpha);/* nb-1 <= |alpha| < nb */

	bk = (double *) calloc(nb, sizeof(double));

	N2 = *N;// for some reason need D to be int, not long.
// evaluate the matern!

	/*
	xscale = abs(x)*(sqrt(8*param["rough"])/ param["range"])
	result = ( param["variance"]/(gamma(param["rough"])* 2^(param["rough"]-1)  ) ) *
			( xscale^param["rough"] *
				besselK(xscale , param["rough"]) )
*/

	xscale = sqrt(8 * (*rough)) / *range;
	varscale =  log(*variance)  - lgammafn(*rough ) -  (*rough -1)*M_LN2;

	for(D=0; D < N2; D++) {
		thisx = fabs(distance[D])*xscale;
		distance[D] = exp(varscale + *rough * log(thisx) )*
				bessel_k_ex(thisx, alpha, 1.0, bk);


	if(isnan(distance[D])) // assume distance is very small
		distance[D] = *variance;
	}

	*range = xscale;
	*rough=varscale;

    free(bk);

}


