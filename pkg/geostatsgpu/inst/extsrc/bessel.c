// this isn't used 
// was only written for debugging

#include <R_ext/Print.h>
#include <math.h>

#define GSL_DBL_EPSILON 2.2204460492503131e-16
#define GSL_SUCCESS 0
#define GSL_SQRT_DBL_MAX 1.3407807929942596e+154
#define GSL_EMAXITER 11

struct gsl_sf_result_e10_struct {
  double val;
  double err;
  int    e10;
};
typedef struct gsl_sf_result_e10_struct gsl_sf_result_e10;


void Rtemme_gamma(double *nu, double * g_1pnu, double * g_1mnu, double *g1, double *g2);

int gsl_sf_temme_gamma(const double nu, double * g_1pnu, double * g_1mnu, double * g1, double * g2);

int gsl_sf_bessel_Knu_scaled_e10_e(const double nu, const double x, gsl_sf_result_e10 * result);

void maternConstants(
	const double *param, 
	// shape range variance nugget anisoRatio anisoAngleRadians anisoAngleDegrees
	int *nuround,
	double *mu,
	double *g_1pnu, 
	double *g_1mnu, 
	double *g1, 
	double *g2,
	double *sinrat,
	double *sinAngle,
	double *cosAngle,
	double *anisoRatioSq,
	double *varscale,
	double *logxscale,
	double *totalVariance
	) {

	*nuround = round(param[0]+0.5);
	*mu = param[0] - *nuround;

	double pi_nu = M_PI * *mu;	
	*sinrat = (fabs(pi_nu) < GSL_DBL_EPSILON ? 1.0 : pi_nu/sin(pi_nu));

	Rtemme_gamma(mu, g_1pnu, g_1mnu, g1, g2);

	*totalVariance = param[3] + param[2];
	*anisoRatioSq = (param[4])*(param[4]);

	*varscale = log(param[2])  - Rf_lgammafn(param[0]) -  (param[0]-1)*M_LN2;
	*logxscale = 1.5 * M_LN2 + 0.5 * log(param[0]) - log(param[1]);
}

void bessel_Knu_scaled(double *nu, double *x, double *result, int *N) {
	int D;
	gsl_sf_result_e10 resG;
	for(D=0;D<*N;++D){

		gsl_sf_bessel_Knu_scaled_e10_e(*nu, x[D], &resG);

		result[D] = (resG.val);
	}
}


int
gsl_sf_bessel_K_scaled_temme(const double mu, const double x,
		double * K_nu, double * K_nup1, double * Kp_nu);


int gsl_sf_bessel_K_scaled_temmeP(double nu, const double x,
                             double * K_nu, double * K_nup1, double * Kp_nu)
{
  int max_iter = 15000;

  const double half_x    = 0.5 * x;
  const double ln_half_x = log(half_x);
  const double half_x_nu = exp(nu*ln_half_x);
  const double pi_nu   = M_PI * nu;
  const double sigma   = -nu * ln_half_x;
  const double sinrat  = (fabs(pi_nu) < GSL_DBL_EPSILON ? 1.0 : pi_nu/sin(pi_nu));
  const double sinhrat = (fabs(sigma) < GSL_DBL_EPSILON ? 1.0 : sinh(sigma)/sigma);
  const double ex = exp(x);

  double sum0, sum1;
  double fk, pk, qk, hk, ck;
  int k = 0;
  int stat_iter;

  double g_1pnu, g_1mnu, g1, g2;
  int stat_g = 0;
  Rtemme_gamma(&nu, &g_1pnu, &g_1mnu, &g1, &g2);

//  Rprintf("bb %d %f %f %f %f %f\n", -1, nu,  g_1pnu, g_1mnu, g1, g2);
  fk = sinrat * (cosh(sigma)*g1 - sinhrat*ln_half_x*g2);
//  Rprintf("cc0 %d %f %f %f %f\n", -1, sinrat, sigma, g1, g2);

  pk = 0.5/half_x_nu * g_1pnu;
  qk = 0.5*half_x_nu * g_1mnu;
  hk = pk;
  ck = 1.0;
  sum0 = fk;
  sum1 = hk;
//  Rprintf("cc1 %d %f %f %f %f\n", -1, sigma, sinhrat, fk, pk);
  while(k < max_iter) {
    double del0;
    double del1;
    k++;
    fk  = (k*fk + pk + qk)/(k*k-nu*nu);
    ck *= half_x*half_x/k;
//    Rprintf("cc2 %d %f %f %f %f\n", -1, fk, pk, ck, qk);
    pk /= (k - nu);
    qk /= (k + nu);
    hk  = -k*fk + pk;
    del0 = ck * fk;
    del1 = ck * hk;
    sum0 += del0;
    sum1 += del1;
    if(fabs(del0) < 0.5*fabs(sum0)*GSL_DBL_EPSILON) break;
  }
//  Rprintf("dd %d %f %f %f %f\n", -1, sum0, sum1, ck, fk);

  *K_nu   = sum0 * ex;
  *K_nup1 = sum1 * 2.0/x * ex;
  *Kp_nu  = - *K_nup1 + nu/x * *K_nu;

  stat_iter = ( k == max_iter ? GSL_EMAXITER : GSL_SUCCESS );
  return 1;//GSL_ERROR_SELECT_2(stat_iter, stat_g);
}


void bessel_K_temme(int *nuround, double *mu, double *x, double *result, int *Nx) {

	int D, Dn, e10;
	double K_mu, K_mup1, Kp_mu;
	double K_nu, K_nup1, K_num1;

	for(D=0;D<*Nx;++D) {

		gsl_sf_bessel_K_scaled_temmeP(*mu, x[D], &K_mu, &K_mup1, &Kp_mu);
//		Rprintf("a %d %f %f %f %f\n", D, x[D], K_mu, K_mup1, Kp_mu);

		K_nu   = K_mu;
		K_nup1 = K_mup1;
		//		for(Dn=0; Dn<*nuround; Dn++) {
		for(Dn=0; Dn<1; Dn++) {
			K_num1 = K_nu;
			K_nu   = K_nup1;
			/* rescale the recurrence to avoid overflow */
#ifdef UNDEF
			if (fabs(K_nu) > GSL_SQRT_DBL_MAX) {
				double p = floor(log(fabs(K_nu))/M_LN10);
				double factor = pow(10.0, p);
				K_num1 /= factor;
				K_nu /= factor;
				e10 += p;
			};
#endif
//			Rprintf("b %f %f %f %f\n", K_num1, K_nu, K_nup1, *mu);
			K_nup1 = 2.0*(*mu+Dn+1)/x[D] * K_nu + K_num1;
			result[D] = K_nup1;
		}
		//		result[D] = K_nu;
	}
}

void bessel_K_p(int *nuround, double *mu, double *nu, double *x, double *result, int *Nx) {

	int D, Dn, e10;
	double K_mu, K_mup1, Kp_mu;
	double K_nu, K_nup1, K_num1, Kp_nu;

	double half_x, ln_half_x, half_x_nu, sigma, sinhrat,ex;

	double sum0, sum1;
	double fk, pk, qk, hk, ck;
	int k;
	int stat_iter;
	double del0;
	double del1;

	double g_1pnu, g_1mnu, g1, g2;
	const double pi_nu = M_PI * (*mu);
	double sinrat =(fabs(pi_nu) < GSL_DBL_EPSILON ? 1.0 : pi_nu/sin(pi_nu));

	int max_iter = 15000;
	Rtemme_gamma(mu, &g_1pnu, &g_1mnu, &g1, &g2);

	for(D=0;D<*Nx;++D) {
		half_x    = 0.5 * x[D];
		ln_half_x = log(half_x);

		half_x_nu = exp((*mu)*ln_half_x);
		sigma   = - (*mu) * ln_half_x;
		sinhrat = (fabs(sigma) < GSL_DBL_EPSILON ? 1.0 : sinh(sigma)/sigma);
//	    Rprintf("b %d %f %f %f %f %f %f\n", D, *mu, *nu,  g_1pnu, g_1mnu, g1, g2);
		fk = sinrat * (cosh(sigma)*g1 - sinhrat*ln_half_x*g2);
//		Rprintf("c0 %d %f %f %f %f\n", D, sinrat, sigma, g1, g2);

		// no need for sigma anymore
		pk = 0.5/half_x_nu * g_1pnu;
		qk = 0.5*half_x_nu * g_1mnu;
		hk = pk;
		ck = 1.0;
		sum0 = fk;
		sum1 = hk;

		// fk, ck, pk, qk hk are vectors
		// keep sum0, sum1, fk, pk, qk, generate kh, ck
		// but qk = (g_1mnu / g_1pnu) * 0.25/pk
		// keep sum0, del0, sum1, fk, log_half_x, x? generate pk, qk, ck?
		del0=0.0;
		del1=0.0;
		k = 0;
		while(k < max_iter) {
			k++;
			fk  = (k*fk + pk + qk)/(k*k-(*mu)* (*mu));
			ck *= half_x*half_x/k; // log(ck) = k * log(x/2) - sumk
			del0 = ck * fk;
			sum0 += del0;

			pk /= (k - (*mu));
			qk /= (k + (*mu));

			hk  = -k*fk + pk;
			del1 = ck * hk; // = ck * (pk - k fk)
			sum1 += del1;
			if(fabs(del0) < 0.5*fabs(sum0)*GSL_DBL_EPSILON) break;
		}
//	  Rprintf("d %d %f %f %f %f\n", D, sum0, pk, ck, fk);

		// Add the matern stuff (on log scale) here?
		ex = exp(x[D]);
		K_mu   = sum0 * ex;
		K_mup1 = sum1 * 2.0/x[D] * ex;
		// is Kp_mu needed?
		Kp_mu  = - K_mup1 + *nu/x[D] * K_mu;

		stat_iter = ( k == max_iter ? GSL_EMAXITER : GSL_SUCCESS );

//		Rprintf("a %d %f %f %f %f\n", D, x[D], K_mu, K_mup1, Kp_mu);

		/* recurse forward to obtain K_num1, K_nu */
		K_nu   = K_mu;
		K_nup1 = K_mup1;

		for(Dn=0; Dn<*nuround; Dn++) {
			K_num1 = K_nu;
			K_nu   = K_nup1;
			/* rescale the recurrence to avoid overflow
			 * not needed if we've scaled to a matern?
			 * */
#ifdef UNDEF
			if (fabs(K_nu) > GSL_SQRT_DBL_MAX) {
				double p = floor(log(fabs(K_nu))/M_LN10);
				Rprintf("%f %f %f\n", K_num1, K_nu, p);
				double factor = pow(10.0, p);
				K_num1 /= factor;
				K_nu /= factor;
				e10 += p;
			};
#endif
			// does this need modifying if we're doing the matern?
			K_nup1 = 2.0*(*mu+Dn+1)/x[D] * K_nu + K_num1;
		}
//		Rprintf("bb %d %f %f %f\n", D, K_nu, K_nup1, *mu);

		result[D] = K_nu;
	}
}
