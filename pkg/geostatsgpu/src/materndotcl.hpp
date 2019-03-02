
std::string maternCLstring = "__kernel void maternCLD("
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
	"__global double *result) {"

	// Get the index of the elements to be processed
	"int Dcell, Drow, Dcol, k;"

	"double dist[2], distRotate[2], distSq;"
	"double logthisx, K_mu, K_mup1, logck, K_nu, K_nup1, Kp_nu, K_num1;"
	"double twoLnHalfX, sigma, half_x_nu, sinhrat;"
	"double ln_half_x, thisx, maternBit, expMaternBit;"

	"double bi, delhi, hi, di, qi, qip1, ai, a1, ci, Qi, s, dels, tmp;"
	"double fk, pk, qk, hk, sum0, sum1, del0, ck, del1;"

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

	"for(Dcell = get_global_id(0); Dcell < Ncell; Dcell += get_global_size(0)) {"

	"	Drow = ceil(0.5 + sqrt(0.25 + 2.0*(Dcell+1) ) ) - 1;"
	"	Dcol = Dcell - round(Drow * (Drow - 1.0) / 2.0);"


	"	dist[0] = coords[Drow*sizeCoordsPadCol] - coords[Dcol*sizeCoordsPadCol];"
	"	dist[1] = coords[Dcol*sizeCoordsPadCol +1] - coords[Drow*sizeCoordsPadCol +1];"
	"	distRotate[0] = costheta *dist[0] + sintheta * dist[1];"
	"	distRotate[1] = sintheta *dist[0] - costheta * dist[1];"
	"	distSq = distRotate[0]*distRotate[0] + distRotate[1]*distRotate[1]/anisoRatioSq;"

	"	logthisx = log(distSq)/2 + logxscale;"

	"	ln_half_x = logthisx - M_LN2;"
	"	thisx = exp(logthisx);"
	"	maternBit = varscale + nu * logthisx;"
	"	expMaternBit = exp(maternBit);"

		// gsl_sf_bessel_K_scaled_temme x < 2
		// gsl_sf_bessel_K_scaled_steed_temme_CF2 x > 2
//		result[Dcol * sizeResultPadRow + Drow] = logthisx;
	"	if(logthisx > 2.0) {" // gsl_sf_bessel_K_scaled_steed_temme_CF2
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

	"		s = 1.0 + Qi*delhi;"

	"		for(k=2; k<=maxIter; k++) {"
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
	"			}" // k loop

	"			hi *= -a1;"
	"		K_nu = exp(-thisx) * exp(maternBit) * sqrt(M_PI/(2.0*thisx)) / s;"//  sqrt(pi)/2 sqrt(2/x)/s =

	"		K_nup1 = K_nu * (mu + thisx + 0.5 - hi)/thisx;"
	"		Kp_nu  = - K_nup1 + mu/thisx * K_nu;"

	"		} else { "// if short distance gsl_sf_bessel_K_scaled_temme
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
	"			del0 = fabs(sum0)+100;"

	"			while( (k < maxIter) && ( fabs(del0) > (epsilon * fabs(sum0)) ) ) {"
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
	"			}" //while loop
			//		result[Dcol * sizeResultPadRow + Drow] = -k;
	"			K_nu   = sum0;"
	"			K_nup1 = sum1 * exp( - ln_half_x);"
	"		}" // short distance

	"		for(k=0; k<nuround; k++) {"
	"			K_num1 = K_nu;"
	"			K_nu   = K_nup1;"
	"			K_nup1 = exp(log(mup1+k) - ln_half_x) * K_nu + K_num1;"
	"		}"
	"	result[Dcol + Drow * sizeResultPadRow] = K_nu;" // lower triangle
	"	result[Drow + Dcol * sizeResultPadRow] = K_nu;" // upper triangle
	"	}" // loop through cells
"};\n\n"//function

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
	"__global float *result) {"
	// Get the index of the elements to be processed
	"int Dcell, Drow, Dcol, k;"

	"float dist[2], distRotate[2], distSq;"
	"float logthisx, K_mu, K_mup1, logck, K_nu, K_nup1, Kp_nu, K_num1;"
	"float twoLnHalfX, sigma, half_x_nu, sinhrat;"
	"float ln_half_x, thisx, maternBit, expMaternBit;"

	"float bi, delhi, hi, di, qi, qip1, ai, a1, ci, Qi, s, dels, tmp;"
	"float fk, pk, qk, hk, sum0, sum1, del0, ck, del1;"

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

	"for(Dcell = get_global_id(0); Dcell < Ncell; Dcell += get_global_size(0)) {"

	"	Drow = ceil(0.5 + sqrt(0.25 + 2.0*(Dcell+1) ) ) - 1;"
	"	Dcol = Dcell - round(Drow * (Drow - 1.0) / 2.0);"


	"	dist[0] = coords[Drow*sizeCoordsPadCol] - coords[Dcol*sizeCoordsPadCol];"
	"	dist[1] = coords[Dcol*sizeCoordsPadCol +1] - coords[Drow*sizeCoordsPadCol +1];"
	"	distRotate[0] = costheta *dist[0] + sintheta * dist[1];"
	"	distRotate[1] = sintheta *dist[0] - costheta * dist[1];"
	"	distSq = distRotate[0]*distRotate[0] + distRotate[1]*distRotate[1]/anisoRatioSq;"

	"	logthisx = log(distSq)/2 + logxscale;"

	"	ln_half_x = logthisx - M_LN2;"
	"	thisx = exp(logthisx);"
	"	maternBit = varscale + nu * logthisx;"
	"	expMaternBit = exp(maternBit);"

		// gsl_sf_bessel_K_scaled_temme x < 2
		// gsl_sf_bessel_K_scaled_steed_temme_CF2 x > 2
//		result[Dcol * sizeResultPadRow + Drow] = logthisx;
	"	if(logthisx > 2.0) {" // gsl_sf_bessel_K_scaled_steed_temme_CF2
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

	"		s = 1.0 + Qi*delhi;"

	"		for(k=2; k<=maxIter; k++) {"
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
	"			}" // k loop

	"			hi *= -a1;"
	"		K_nu = exp(-thisx) * exp(maternBit) * sqrt(M_PI/(2.0*thisx)) / s;"//  sqrt(pi)/2 sqrt(2/x)/s =

	"		K_nup1 = K_nu * (mu + thisx + 0.5 - hi)/thisx;"
	"		Kp_nu  = - K_nup1 + mu/thisx * K_nu;"

	"		} else { "// if short distance gsl_sf_bessel_K_scaled_temme
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
	"			del0 = fabs(sum0)+100;"

	"			while( (k < maxIter) && ( fabs(del0) > (epsilon * fabs(sum0)) ) ) {"
	"				k++;"
	"				logck += twoLnHalfX - log((float)k);"
	"				ck = exp(logck);"
	"				fk  = (k*fk + pk + qk)/(k*k-muSq);"

	"				del0 = ck * fk;"
	"				sum0 += del0;"

	"				pk /= (k - (mu));"
	"				qk /= (k + (mu));"

	"				hk  = -k*fk + pk;"
	"				del1 = ck * hk; "	
	"				sum1 += del1;"
	"			}" //while loop
			//		result[Dcol * sizeResultPadRow + Drow] = -k;
	"			K_nu   = sum0;"
	"			K_nup1 = sum1 * exp( - ln_half_x);"
	"		}" // short distance

	"		for(k=0; k<nuround; k++) {"
	"			K_num1 = K_nu;"
	"			K_nu   = K_nup1;"
	"			K_nup1 = exp(log(mup1+k) - ln_half_x) * K_nu + K_num1;"
	"		}"
	"	result[Dcol + Drow * sizeResultPadRow] = K_nu;" // lower triangle
	"	result[Drow + Dcol * sizeResultPadRow] = K_nu;" // upper triangle
	"	}" // loop through cells
"};\n"
;


