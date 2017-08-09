
x = seq(1, 3, len=2)
nu = 1.5
nuround = round(nu+0.5)
mu = nu - nuround

y1 = besselK(x, nu)
y2 = besselK(x, nu, TRUE)
range(log(y2) - x - log(y1), na.rm=T)

x.vec = x
nu.vec = rep(nu, length(x))
jj <- .C("bessel_Knu_scaled_e", as.double(nu.vec), as.double(x.vec), 
    as.integer(length(x.vec)), val = as.double(x.vec), err = as.double(x.vec), 
    status = as.integer(0 * x.vec), PACKAGE = "gsl")
y3 = jj$val

range(y2-y3, na.rm=T)


y4 = .C("bessel_Knu_scaled", nu, as.double(x), as.double(x), 
    as.integer(length(x)), PACKAGE='geostatsgpu')[[3]] 

range(y3-y4, na.rm=TRUE)

y5 = .C("bessel_K_temme", as.integer(nuround), as.double(mu), as.double(x), as.double(x), 
    as.integer(length(x)), PACKAGE='geostatsgpu')[[4]] 


range(y5-y4, na.rm=TRUE)

y6 = .C("bessel_K_p", 
    as.integer(nuround), as.double(nu), as.double(mu), 
    as.double(x), as.double(x), 
    as.integer(length(x)), PACKAGE='geostatsgpu')[[5]] 

range(y5-y6, na.rm=TRUE)




.C("Rtemme_gamma", 
    as.double(nu), 
    as.double(0), as.double(0), 
    as.double(0), as.double(0), PACKAGE = "geostatsgpu")



#int stat_g = gsl_sf_temme_gamma(nu, &g_1pnu, &g_1mnu, &g1, &g2);
#gsl_sf_temme_gamma(const double nu, double * g_1pnu, double * g_1mnu, double * g1, double * g2)
#{
anu = abs(nu)#;    /* functions are even */
xHere = 4.0*anu - 1.0;

  gsl_sf_result r_g1;
  gsl_sf_result r_g2;
  cheb_eval_e(&g1_cs, x, &r_g1);
  cheb_eval_e(&g2_cs, x, &r_g2);
  *g1 = r_g1.val;
  *g2 = r_g2.val;
  *g_1mnu = 1.0/(r_g2.val + nu * r_g1.val);
  *g_1pnu = 1.0/(r_g2.val - nu * r_g1.val);
#  return GSL_SUCCESS;
#}


sum0 = sum1 = 1.0
fk = pk = qk = hk = ck= 1.0
k = 0L
stat_iter = 0L

max_iter = 15000;
pi_nu   = pi * nu;
sinrat  = pi_nu/sin(pi_nu);

ex = exp(x);
half_x = 0.5 * x;
ln_half_x = log(half_x);
half_x_nu = exp(nu*ln_half_x);
sigma   = -nu * ln_half_x;
sinhrat = sinh(sigma)/sigma;




fk = sinrat * (cosh(sigma)*g1 - sinhrat*ln_half_x*g2);
pk = 0.5/half_x_nu * g_1pnu;
qk = 0.5*half_x_nu * g_1mnu;
hk = pk;
ck = 1.0;
sum0 = fk;
sum1 = hk;
while(k < max_iter) {
  double del0;
  double del1;
  k++;
  fk  = (k*fk + pk + qk)/(k*k-nu*nu);
  ck *= half_x*half_x/k;
  pk /= (k - nu);
  qk /= (k + nu);
  hk  = -k*fk + pk;
  del0 = ck * fk;
  del1 = ck * hk;
  sum0 += del0;
  sum1 += del1;
  if(fabs(del0) < 0.5*fabs(sum0)*GSL_DBL_EPSILON) break;
}

*K_nu   = sum0 * ex;
*K_nup1 = sum1 * 2.0/x * ex;
*Kp_nu  = - *K_nup1 + nu/x * *K_nu;

stat_iter = ( k == max_iter ? GSL_EMAXITER : GSL_SUCCESS );
return GSL_ERROR_SELECT_2(stat_iter, stat_g);
}