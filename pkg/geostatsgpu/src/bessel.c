#include <R_ext/Print.h>
#include <math.h>
#include <stdlib.h>

#define GSL_DBL_EPSILON 2.2204460492503131e-16
#define GSL_SUCCESS 0
#define GSL_SQRT_DBL_MAX 1.3407807929942596e+154
#define GSL_EMAXITER 11

struct gsl_sf_result_struct {
  double val;
  double err;
};
typedef struct gsl_sf_result_struct gsl_sf_result;


struct gsl_sf_result_e10_struct {
  double val;
  double err;
  int    e10;
};
typedef struct gsl_sf_result_e10_struct gsl_sf_result_e10;

struct gsl_cheb_series_struct {

  double * c;   /* coefficients                */
  size_t order; /* order of expansion          */
  double a;     /* lower interval point        */
  double b;     /* upper interval point        */

  size_t order_sp;

  double * f;   /* function evaluated at chebyschev points  */
};
typedef struct gsl_cheb_series_struct gsl_cheb_series;



/* nu = (x+1)/4, -1<x<1, 1/(2nu)(1/Gamma[1-nu]-1/Gamma[1+nu]) */
static double g1_dat[14] = {
		-1.14516408366268311786898152867,
		0.00636085311347084238122955495,
		0.00186245193007206848934643657,
		0.000152833085873453507081227824,
		0.000017017464011802038795324732,
		-6.4597502923347254354668326451e-07,
		-5.1819848432519380894104312968e-08,
		4.5189092894858183051123180797e-10,
		3.2433227371020873043666259180e-11,
		6.8309434024947522875432400828e-13,
		2.8353502755172101513119628130e-14,
		-7.9883905769323592875638087541e-16,
		-3.3726677300771949833341213457e-17,
		-3.6586334809210520744054437104e-20
};

static gsl_cheb_series g1_cs = {
		g1_dat,
		13,
		-1, 1,
		7
};

/* nu = (x+1)/4, -1<x<1,  1/2 (1/Gamma[1-nu]+1/Gamma[1+nu]) */
static double g2_dat[15] =
{
		1.882645524949671835019616975350,
		-0.077490658396167518329547945212,
		-0.018256714847324929419579340950,
		0.0006338030209074895795923971731,
		0.0000762290543508729021194461175,
		-9.5501647561720443519853993526e-07,
		-8.8927268107886351912431512955e-08,
		-1.9521334772319613740511880132e-09,
		-9.4003052735885162111769579771e-11,
		4.6875133849532393179290879101e-12,
		2.2658535746925759582447545145e-13,
		-1.1725509698488015111878735251e-15,
		-7.0441338200245222530843155877e-17,
		-2.4377878310107693650659740228e-18,
		-7.5225243218253901727164675011e-20
};


static gsl_cheb_series g2_cs = {
		g2_dat,
		14,
		-1, 1,
		8
};


int gsl_sf_temme_gamma(const double nu, double * g_1pnu, double * g_1mnu, double * g1, double * g2);

int gsl_sf_bessel_Knu_scaled_e10_e(const double nu, const double x, gsl_sf_result_e10 * result);

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





static inline int
cheb_eval_e(const gsl_cheb_series * cs,
		const double x,
		gsl_sf_result * result)
{
	int j;
	double d  = 0.0;
	double dd = 0.0;

	double y  = (2.0*x - cs->a - cs->b) / (cs->b - cs->a);
	double y2 = 2.0 * y;

	double e = 0.0;

	for(j = cs->order; j>=1; j--) {
		double temp = d;
		d = y2*d - dd + cs->c[j];
		e += fabs(y2*temp) + fabs(dd) + fabs(cs->c[j]);
		dd = temp;
	}

	{
		double temp = d;
		d = y*d - dd + 0.5 * cs->c[0];
		e += fabs(y*temp) + fabs(dd) + 0.5 * fabs(cs->c[0]);
	}

	result->val = d;
	result->err = GSL_DBL_EPSILON * e + fabs(cs->c[cs->order]);

	return GSL_SUCCESS;
}

void Rtemme_gamma(double *nu, double * g_1pnu, double * g_1mnu, double *g1, double *g2) {

	const double anu = fabs(*nu);    /* functions are even */
	const double x = 4.0*anu - 1.0;
	gsl_sf_result r_g1;
	gsl_sf_result r_g2;

	cheb_eval_e(&g1_cs, x, &r_g1);
	cheb_eval_e(&g2_cs, x, &r_g2);

	*g1 = r_g1.val;
	*g2 = r_g2.val;
	*g_1mnu = 1.0/(r_g2.val + *nu * r_g1.val);
	*g_1pnu = 1.0/(r_g2.val - *nu * r_g1.val);
}

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
