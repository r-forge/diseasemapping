// TO DO 
// create float string using maternClcommonString
// reduce the number of variables so they all fit in private memory
// use macros to replace different variables in the if(logthisx > 2.0) blocks
// #define sigma hi and get rid of sigma as a variable?
// make globalSize local?

std::string maternCLcommonString =
  
// matrix stacked row-wise, excluding diagonal
// i = row, j  = col, k = cell
// k = i^2/2 - i/2 - i + j + 1
// if j = i-1
// i^2/2 - 1/2 i -k = 0 or i = [1/2 pm sqrt (1/4 + 4 * (1/2) * k)] / 2 * (1/2)  
//   i = 1/2 + sqrt (1/4 + 2*k)
// kk = 1:12;ii = ceiling(1/2 + sqrt(0.25 + 2*kk));ii
// jj = kk - 1 - ii^2/2 + (3/2)*ii;jj
// Drow = ii-1;jj = kk - 1 - Drow-1^2/2 + (3/2)*ii;jj
// Dcell +1 = (Drow+1)*Drow/2 - Drow + Dcol + 1 
// Dcell = Drow (Drow/2 - 1/2) + Dcol
// Dcol = Dcell - Drow * (Drow-1)/2
// Ncell = (N-1) * N / 2


"const int globalSize = get_global_size(0);"
"int Dcell, Drow, Dcol, k;\n"
"for(Dcell = get_global_id(0); Dcell < Ncell; Dcell += globalSize) {\n"

"	Drow = ceil(0.5 + sqrt(0.25 + 2.0*(Dcell+1) ) ) - 1;\n"
"	Dcol = Dcell - round(Drow * (Drow - 1.0) / 2.0);\n"


"	dist[0] = coords[Drow*sizeCoordsPadCol] - coords[Dcol*sizeCoordsPadCol];"
"	dist[1] = coords[Dcol*sizeCoordsPadCol +1] - coords[Drow*sizeCoordsPadCol +1];"
"	distRotate[0] = costheta *dist[0] + sintheta * dist[1];"
"	distRotate[1] = sintheta *dist[0] - costheta * dist[1];"
"	distSq = distRotate[0]*distRotate[0] + distRotate[1]*distRotate[1]/anisoRatioSq;"

"	logthisx = log(distSq)/2 + logxscale;"

"	ln_half_x = logthisx - logTwo;"
"	thisx = exp(logthisx);"
"	maternBit = varscale + nu * logthisx;"
"	expMaternBit = exp(maternBit);\n"
// gsl_sf_bessel_K_scaled_temme x < 2
// gsl_sf_bessel_K_scaled_steed_temme_CF2 x > 2
//		result[Dcol * sizeResultPadRow + Drow] = logthisx;
"	if(logthisx > 2.0) {\n" // gsl_sf_bessel_K_scaled_steed_temme_CF2
"		K_nu = 0.0;"
"		bi = 2.0*(1.0 + thisx);"
"		di = 1.0/bi;"

"		delhi = di;"
"		hi = delhi;"// was di, now di/exp(maternBit)

"		qi   = 0.0;"
"		qip1 = 1.0;"

"		ai = -(0.25 - muSq);"
"		a1 = ai;"
"		ci = -ai;"
"		Qi = -ai;"

"		s = 1.0 + Qi*delhi;\n"

"		for(k=2; k<=maxIter; k++) {\n"
"			ai -= 2.0*(k-1);"
"			ci  = -ai*ci/k;"
"			tmp  = (qi - bi*qip1)/ai;"
"			qi   = qip1;"
"			qip1 = tmp;"
"			Qi += ci*qip1;"
"			bi += 2.0;"
"			di  = 1.0/(bi + ai*di);"
"			delhi = (bi*di - 1.0) * delhi;"
"			hi += delhi;"
"			dels = Qi*delhi;"
"			s += dels;"
"			if(fabs(dels/s) < epsilon) break;"
"			}\n" // k loop
"			hi *= -a1;"
//"		K_nu = exp(-thisx) * exp(maternBit) * sqrt(M_PI/(2.0*thisx)) / s;"//  sqrt(pi)/2 sqrt(2/x)/s =
"		K_nu = exp(logSqrtHalfPi + maternBit - thisx - logthisx / 2) / s;"//  sqrt(pi)/2 sqrt(2/x)/s =

"		K_nup1 = K_nu * (mu + thisx + 0.5 - hi)/thisx;"
"		Kp_nu  = - K_nup1 + mu/thisx * K_nu;\n"

"		} else {\n"// if short distance gsl_sf_bessel_K_scaled_temme
"			twoLnHalfX = 2*ln_half_x;"
"			sigma   = - mu * ln_half_x;"
"			half_x_nu = exp(-sigma);"
"			sinhrat = sinh(sigma)/sigma;"

"			fk = sinrat * (cosh(sigma)*g1 - sinhrat*ln_half_x*g2);"
"			pk = 0.5/half_x_nu * g_1pnu;"
"			qk = 0.5*half_x_nu * g_1mnu;"
"			hk = pk;"
"			sum0 = fk*expMaternBit;"
"			sum1 = hk*expMaternBit;"
"			k=0;"
"			logck = maternBit;"
"			del0 = fabs(sum0)+100;\n"
"			while( (k < maxIter) && ( fabs(del0) > (epsilon * fabs(sum0)) ) ) {\n"
"				k++;"
"				logck += twoLnHalfX - log((double)k);"
"				ck = exp(logck);"
"				fk  = (k*fk + pk + qk)/(k*k-muSq);"

"				del0 = ck * fk;"
"				sum0 += del0;"

"				pk /= (k - (mu));"
"				qk /= (k + (mu));"

"				hk  = -k*fk + pk;"
"				del1 = ck * hk; "	
"				sum1 += del1;"
"			}\n" //while loop
//		result[Dcol * sizeResultPadRow + Drow] = -k;
"			K_nu   = sum0;"
"			K_nup1 = sum1 * exp( - ln_half_x);"
"		}\n" // short distance
"		for(k=0; k<nuround; k++) {\n"
"			K_num1 = K_nu;"
"			K_nu   = K_nup1;"
"			K_nup1 = exp(log(mup1+k) - ln_half_x) * K_nu + K_num1;"
"		}\n"
"#ifdef assignLower\n"
"	  result[Dcol + Drow * sizeResultPadRow] = K_nu;\n" // lower triangle
"# endif\n"
"#ifdef assignUpper\n;"
"	result[Drow + Dcol * sizeResultPadRow] = K_nu;\n"//K_nu;\n" // upper triangle
"# endif\n"
"	}\n"; // loop through cells



// note that either assignLower or assignUpper or both must be defined
// otherwise the values won't be saved
std::string maternCLstring =
  "\n#pragma OPENCL EXTENSION cl_khr_fp64: enable\n\n"
  "\n#define logSqrtHalfPi 0.2257913526447273278031\n"
  "\n#define logTwo M_LN2\n"
  
  "\n__kernel void maternCLD("
  "const unsigned int Ncell,"
  "const unsigned int sizeCoordsPadCol,"
  "const unsigned int sizeResultPadRow,"
  "const unsigned int sizeResultPadCol,"
  "const unsigned int maxIter,"
  "const double nu,"
  "const int nuround,"
  "const double mu,"
  "const double muSq,"
  "const double mup1,"
  "const double costheta,"
  "const double sintheta,"
  "const double anisoRatioSq,"
  "const double varscale,"
  "const double logxscale,"
  "const double sinrat,"
  "const double g_1pnu,"
  "const double g_1mnu,"
  "const double g1,"
  "const double g2,"
  "const double epsilon,"
  "__global const double *coords,"
  "__global double *result) {\n"
  
  "double dist[2], distRotate[2], distSq,"
  "  logthisx, K_mu, K_mup1, logck, K_nu, K_nup1, Kp_nu, K_num1,"
  "  twoLnHalfX, sigma, half_x_nu, sinhrat,"
  "  ln_half_x, thisx, maternBit, expMaternBit,"
  "  bi, delhi, hi, di, qi, qip1, ai, a1, ci, Qi, s, dels, tmp,"
  "  fk, pk, qk, hk, sum0, sum1, del0, ck, del1;\n"  
  + maternCLcommonString +
  "};\n\n";//function


std::string maternCLstringF = 
"__kernel void maternCLF("
	"const unsigned int Ncell,"
	"const unsigned int sizeCoordsPadCol,"
	"const unsigned int sizeResultPadRow,"
	"const unsigned int sizeResultPadCol,"
	"const unsigned int maxIter,"
	"const float nu,"
	"const int nuround,"
	"const float mu,"
	"const float muSq,"
	"const float mup1,"
	"const float costheta,"
	"const float sintheta,"
	"const float anisoRatioSq,"
	"const float varscale,"
	"const float logxscale,"
	"const float sinrat,"
	"const float g_1pnu,"
	"const float g_1mnu,"
	"const float g1,"
	"const float g2,"
	"const float epsilon,"
	"__global float *coords,"
	"__global float *result) {\n"
	// Get the index of the elements to be processed
	"int Dcell, Drow, Dcol, k;\n"

	"float dist[2], distRotate[2], distSq;\n"
	"float logthisx, K_mu, K_mup1, logck, K_nu, K_nup1, Kp_nu, K_num1;\n"
	"float twoLnHalfX, sigma, half_x_nu, sinhrat;\n"
	"float ln_half_x, thisx, maternBit, expMaternBit;\n"

	"float bi, delhi, hi, di, qi, qip1, ai, a1, ci, Qi, s, dels, tmp;\n"
	"float fk, pk, qk, hk, sum0, sum1, del0, ck, del1;\n"

	// matrix stacked row-wise, excluding diagonal
	// i = row, j  = col, k = cell
	// k = i^2/2 - i/2 - i + j + 1
	// if j = i-1
	// i^2/2 - 1/2 i -k = 0 or i = [1/2 pm sqrt (1/4 + 4 * (1/2) * k)] / 2 * (1/2)  
	//   i = 1/2 + sqrt (1/4 + 2*k)
	// kk = 1:12;ii = ceiling(1/2 + sqrt(0.25 + 2*kk));ii
	// jj = kk - 1 - ii^2/2 + (3/2)*ii;jj
	// Drow = ii-1;jj = kk - 1 - Drow-1^2/2 + (3/2)*ii;jj
	// Dcell +1 = (Drow+1)*Drow/2 - Drow + Dcol + 1 
	// Dcell = Drow (Drow/2 - 1/2) + Dcol
	// Dcol = Dcell - Drow * (Drow-1)/2
	// Ncell = (N-1) * N / 2

	"for(Dcell = get_global_id(0); Dcell < Ncell; Dcell += get_global_size(0)) {\n"

	"	Drow = ceil(0.5 + sqrt(0.25 + 2.0*(Dcell+1) ) ) - 1;\n"
	"	Dcol = Dcell - round(Drow * (Drow - 1.0) / 2.0);\n"


	"	dist[0] = coords[Drow*sizeCoordsPadCol] - coords[Dcol*sizeCoordsPadCol];\n"
	"	dist[1] = coords[Dcol*sizeCoordsPadCol +1] - coords[Drow*sizeCoordsPadCol +1];\n"
	"	distRotate[0] = costheta *dist[0] + sintheta * dist[1];\n"
	"	distRotate[1] = sintheta *dist[0] - costheta * dist[1];\n"
	"	distSq = distRotate[0]*distRotate[0] + distRotate[1]*distRotate[1]/anisoRatioSq;\n"

	"	logthisx = log(distSq)/2 + logxscale;\n"

	"	ln_half_x = logthisx - M_LN2;\n"
	"	thisx = exp(logthisx);\n"
	"	maternBit = varscale + nu * logthisx;\n"
	"	expMaternBit = exp(maternBit);\n"

		// gsl_sf_bessel_K_scaled_temme x < 2
		// gsl_sf_bessel_K_scaled_steed_temme_CF2 x > 2
//		result[Dcol * sizeResultPadRow + Drow] = logthisx;
	"	if(logthisx > 2.0) {\n" // gsl_sf_bessel_K_scaled_steed_temme_CF2
	"		K_nu = 0.0;\n"
	"		bi = 2.0*(1.0 + thisx);\n"
	"		di = 1.0/bi;\n"

	"		delhi = di;\n"
	"		hi = delhi;\n"// was di, now di/exp(maternBit)

	"		qi   = 0.0;\n"
	"		qip1 = 1.0;\n"

	"		ai = -(0.25 - muSq);\n"
	"		a1 = ai;\n"
	"		ci = -ai;\n"
	"		Qi = -ai;\n"

	"		s = 1.0 + Qi*delhi;\n"

	"		for(k=2; k<=maxIter; k++) {\n"
	"			ai -= 2.0*(k-1);\n"
	"			ci  = -ai*ci/k;\n"
	"			tmp  = (qi - bi*qip1)/ai;\n"
	"			qi   = qip1;\n"
	"			qip1 = tmp;\n"
	"			Qi += ci*qip1;\n"
	"			bi += 2.0;\n"
	"			di  = 1.0/(bi + ai*di);\n"
	"			delhi = (bi*di - 1.0) * delhi;\n"
	"			hi += delhi;\n"
	"			dels = Qi*delhi;\n"
	"			s += dels;\n"
	"			if(fabs(dels/s) < epsilon) break;\n"
	"			}\n" // k loop

	"			hi *= -a1;\n"
	"		K_nu = exp(-thisx) * exp(maternBit) * sqrt(M_PI/(2.0*thisx)) / s;\n"//  sqrt(pi)/2 sqrt(2/x)/s =

	"		K_nup1 = K_nu * (mu + thisx + 0.5 - hi)/thisx;\n"
	"		Kp_nu  = - K_nup1 + mu/thisx * K_nu;\n"

	"		} else { \n"// if short distance gsl_sf_bessel_K_scaled_temme
	"			twoLnHalfX = 2*ln_half_x;\n"
	"			sigma   = - mu * ln_half_x;\n"
	"			half_x_nu = exp(-sigma);\n"
	"			sinhrat = sinh(sigma)/sigma;\n"

	"			fk = sinrat * (cosh(sigma)*g1 - sinhrat*ln_half_x*g2);\n"
	"			pk = 0.5/half_x_nu * g_1pnu;\n"
	"			qk = 0.5*half_x_nu * g_1mnu;\n"
	"			hk = pk;\n"
	"			sum0 = fk*expMaternBit;\n"
	"			sum1 = hk*expMaternBit;\n"
	"			k=0;\n"
	"			logck = maternBit;\n"
	"			del0 = fabs(sum0)+100;\n"

	"			while( (k < maxIter) && ( fabs(del0) > (epsilon * fabs(sum0)) ) ) {\n"
	"				k++;\n"
	"				logck += twoLnHalfX - log((float)k);\n"
	"				ck = exp(logck);\n"
	"				fk  = (k*fk + pk + qk)/(k*k-muSq);\n"

	"				del0 = ck * fk;\n"
	"				sum0 += del0;\n"

	"				pk /= (k - (mu));\n"
	"				qk /= (k + (mu));\n"

	"				hk  = -k*fk + pk;\n"
	"				del1 = ck * hk;\n "	
	"				sum1 += del1;\n"
	"			}\n" //while loop
			//		result[Dcol * sizeResultPadRow + Drow] = -k;
	"			K_nu   = sum0;\n"
	"			K_nup1 = sum1 * exp( - ln_half_x);\n"
	"		}\n" // short distance

	"		for(k=0; k<nuround; k++) {\n"
	"			K_num1 = K_nu;\n"
	"			K_nu   = K_nup1;\n"
	"			K_nup1 = exp(log(mup1+k) - ln_half_x) * K_nu + K_num1;\n"
	"		}\n"
	"	result[Dcol + Drow * sizeResultPadRow] = K_nu;\n" // lower triangle
	"	result[Drow + Dcol * sizeResultPadRow] = K_nu;\n" // upper triangle
	"	}\n" // loop through cells
"};\n"
;


