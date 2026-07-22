/* Matérn correlation via Temme series / Steed continued fraction.
   Ported from gpuRandom OpenCL maternBatch (itself based on GSL bessel_Knu).
   Iterates on the Matérn function rather than calling Bessel K separately. */

#include "geostatsp.h"
#include <math.h>
#include <stdlib.h>

#define MATERN_EPSILON 1e-12
#define MATERN_MAXITER 5000
#define GSL_DBL_EPSILON 2.2204460492503131e-16
#define GSL_SUCCESS 0

/* ---- Temme gamma helpers (from gpuRandom/src/besselUtils.c / GSL) ---- */

struct gsl_sf_result_struct {
  double val;
  double err;
};
typedef struct gsl_sf_result_struct gsl_sf_result;

struct gsl_cheb_series_struct {
  double * c;
  size_t order;
  double a;
  double b;
  size_t order_sp;
  double * f;
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
    7,
    NULL
};

/* nu = (x+1)/4, -1<x<1,  1/2 (1/Gamma[1-nu]+1/Gamma[1+nu]) */
static double g2_dat[15] = {
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
    8,
    NULL
};

static int
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

  for(j = (int)cs->order; j>=1; j--) {
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

static void Rtemme_gamma(double nu, double * g_1pnu, double * g_1mnu,
                         double *g1, double *g2)
{
  const double anu = fabs(nu);
  const double x = 4.0*anu - 1.0;
  gsl_sf_result r_g1;
  gsl_sf_result r_g2;

  cheb_eval_e(&g1_cs, x, &r_g1);
  cheb_eval_e(&g2_cs, x, &r_g2);

  *g1 = r_g1.val;
  *g2 = r_g2.val;
  *g_1mnu = 1.0/(r_g2.val + nu * r_g1.val);
  *g_1pnu = 1.0/(r_g2.val - nu * r_g1.val);
}

/* ---- Short-distance Temme series (returns K_nu, K_nup1 via pointers) ---- */

static void maternShort(
    double ln_half_x,
    double maternBit,
    double expMaternBit,
    double mu,
    double muSq,
    double sinrat,
    double g1,
    double g2,
    double g_1pnu,
    double g_1mnu,
    double *K_nu,
    double *K_nup1)
{
  double twoLnHalfX, sigma, sinhrat, fk, pk, qk, hk, half_x_nu;
  double sum0, sum1, ck, logck, del0, del1;
  int k;

  twoLnHalfX = 2.0 * ln_half_x;
  sigma   = - mu * ln_half_x;
  half_x_nu = exp(-sigma);
  sinhrat = sinh(sigma)/sigma;
  if(ISNAN(sinhrat)) sinhrat = 1.0;

  fk = sinrat * (cosh(sigma)*g1 - sinhrat*ln_half_x*g2);
  pk = 0.5/half_x_nu * g_1pnu;
  qk = 0.5*half_x_nu * g_1mnu;
  hk = pk;
  sum0 = fk*expMaternBit;
  sum1 = hk*expMaternBit;
  logck = maternBit;

  /* Require a few terms before testing convergence: the first correction
     can underflow to exactly zero (e.g. shape=1/2, thisx=3) even though
     the series has not yet left the initial fk*expMaternBit value. */
  for(k = 1; k <= MATERN_MAXITER; k++) {
    logck += twoLnHalfX - log((double)k);
    ck = exp(logck);
    fk  = (k*fk + pk + qk)/(k*k - muSq);

    del0 = ck * fk;
    sum0 += del0;

    pk /= (k - mu);
    qk /= (k + mu);

    hk  = -k*fk + pk;
    del1 = ck * hk;
    sum1 += del1;

    if(k >= 3 && fabs(del0) <= (MATERN_EPSILON * fabs(sum0)) &&
       fabs(del1) <= (MATERN_EPSILON * fabs(sum1))) {
      break;
    }
  }

  *K_nu   = sum0;
  *K_nup1 = sum1 * exp( - ln_half_x);
}

/* ---- Long-distance Steed continued fraction ---- */

static void maternLong(
    double thisx,
    double expMaternBit,
    double mu,
    double muSq,
    double *K_nu,
    double *K_nup1)
{
  double bi, di, delhi, hi, qi, qip1, ai, ci, Qi, s, tmp;
  double dels, a1;
  int k;

  bi = 2.0*(1.0 + thisx);
  di = 1.0/bi;
  delhi = di;
  hi = delhi;
  qi   = 0.0;
  qip1 = 1.0;
  ai = -(0.25 - muSq);
  a1 = ai;
  ci = -ai;
  Qi = -ai;
  s = 1.0 + Qi*delhi;

  for(k=2; k<=MATERN_MAXITER; k++) {
    ai -= 2.0*(k-1);
    ci  = -ai*ci/k;
    tmp  = (qi - bi*qip1)/ai;
    qi   = qip1;
    qip1 = tmp;
    Qi += ci*qip1;
    bi += 2.0;
    di  = 1.0/(bi + ai*di);
    delhi = (bi*di - 1.0) * delhi;
    hi += delhi;
    dels = Qi*delhi;
    s += dels;
    if(fabs(dels/s) < MATERN_EPSILON) break;
  }
  hi *= -a1;

  *K_nu = exp(-thisx) * expMaternBit * sqrt(M_PI_2/thisx) / s;
  *K_nup1 = (*K_nu) * (mu + thisx + 0.5 - hi)/thisx;
}

/* ---- Parameter setup (once per call) ---- */

void maternParamsSet(maternParams *p,
                     double range, double shape, double variance)
{
  double pi_nu;

  p->range = range;
  p->shape = shape;
  p->variance = variance;
  p->varscale = log(variance) - lgammafn(shape) - (shape - 1.0)*M_LN2;
  p->logxscale = 1.5 * M_LN2 + 0.5 * log(shape) - log(range);

  p->nuround = (int)(shape + 0.5);
  p->mu = shape - (double)p->nuround;
  p->muSq = p->mu * p->mu;
  p->mup1 = p->mu + 1.0;

  pi_nu = M_PI * p->mu;
  p->sinrat = (fabs(pi_nu) < GSL_DBL_EPSILON) ? 1.0 : pi_nu/sin(pi_nu);

  Rtemme_gamma(p->mu, &p->g_1pnu, &p->g_1mnu, &p->g1, &p->g2);
}

/* ---- Scalar evaluator: distSq is squared Euclidean (anisotropic) distance ---- */

double maternPoint(double distSq, const maternParams *p)
{
  double logthisx, ln_half_x, thisx, maternBit, expMaternBit;
  double K_nu, K_nup1, K_num1;
  int k;
  const double truncate = p->variance * 1e-06;

  if(!(distSq > 0.0)) {
    /* zero or negative distance */
    return p->variance;
  }

  /* gaussian limit for very large shape */
  if(p->nuround >= 1000) {
    return p->variance * exp(-2.0 * distSq / (p->range * p->range));
  }

  logthisx = 0.5 * log(distSq) + p->logxscale;
  thisx = exp(logthisx);

  /* edge cases mirroring previous matern.c behaviour */
  if(ISNAN(thisx)) {
    if(!R_FINITE(p->logxscale)) {
      /* range is probably zero */
      return (distSq < truncate) ? p->variance : 0.0;
    }
    /* finite range, distance must be zero / underflow */
    return p->variance;
  }

  maternBit = p->varscale + p->shape * logthisx;
  expMaternBit = exp(maternBit);
  ln_half_x = logthisx - M_LN2;

  if(logthisx > 1.5) {
    maternLong(thisx, expMaternBit, p->mu, p->muSq, &K_nu, &K_nup1);
  } else {
    maternShort(ln_half_x, maternBit, expMaternBit,
                p->mu, p->muSq, p->sinrat,
                p->g1, p->g2, p->g_1pnu, p->g_1mnu,
                &K_nu, &K_nup1);
  }

  for(k = 0; k < p->nuround; k++) {
    K_num1 = K_nu;
    K_nu   = K_nup1;
    K_nup1 = exp(log(p->mup1 + k) - ln_half_x) * K_nu + K_num1;
  }

  if(ISNAN(K_nu)) {
    return (thisx < 1.0) ? p->variance : 0.0;
  }

  return K_nu;
}
