

__kernel void maternSingleIndexCL(
	__local const unsigned int Ncell,
	__local const unsigned int sizeCoordsPadCol,
	__local const unsigned int sizeResultPadRow,
	__local const unsigned int sizeResultPadCol,
	__local const unsigned int maxIter,
	__local const double nu,
	__local const int nuround,
	__local const double mu,
	__local const double muSq,
	__local const double mup1,
	__local const double costheta,
	__local const double sintheta,
	__local const double anisoRatioSq,
	__local const double varscale,
	__local const double logxscale,
	__local const double sinrat,
	__local const double g_1pnu,
	__local const double g_1mnu,
	__local const double g1,
	__local const double g2,
	__local const double epsilon,
	__local const double *coords,
	__global double *result) {
	// Get the index of the elements to be processed

	unsigned int Dcell, Drow, Dcol, k;

	double dist[2], distRotate[2], distSq;
	double logthisx, K_mu, K_mup1, logck, K_nu, K_nup1, Kp_nu;
	double twoLnHalfX, sigma, half_x_nu, sinhrat;
	double ln_half_x, thisx, maternBit, expMaternBit;

	double bi, delhi, hi, di, qi, qip1, ai, a1, ci, Qi, s, dels, tmp;
	double fk, pk, qk, hk, sum0, sum1, logck, del0, ck, del1;

	// matrix stacked row-wise, excluding diagonal
	// k = i * (i-1) + j
	// i = floor( (1 + sqrt(1 + 8*k)) / 2)
	// j = k - i*(i-1)/2
	// Drow = i-1, Dcol = j-1, Dcell = k-1
	// Ncell = (N-1) * N / 2

	for(Dcell = get_global_id(0), 
		Dcell < Ncell, 
		Dcell += get_global_size(0)) {

		Drow = floor( (1 + sqrt(9 + 8*Dcell)) / 2)-1;
		Dcol = Dcell - (Drow + 1) * Drow / 2;

		dist[0] = coords[Drow*sizeCoordsPadCol] - coords[Dcol*sizeCoordsPadCol];
		dist[1] = coords[Dcol*sizeCoordsPadCol +1] - coords[Drow*sizeCoordsPadCol +1];
		distRotate[0] = costheta *dist[0] - sintheta * dist[1];
		distRotate[1] = sintheta *dist[0] + costheta * dist[1];
		distSq = distRotate[0]*distRotate[0] + distRotate[1]*distRotate[1]/anisoRatioSq;

		//		result[Dcol * sizeResultPadRow + Drow] = distSq; // upper triangle

		logthisx = log(distSq)/2 + logxscale;

		ln_half_x = logthisx - M_LN2;
		thisx = exp(logthisx);
		maternBit = varscale + nu * logthisx;
		expMaternBit = exp(maternBit);

		// gsl_sf_bessel_K_scaled_temme x < 2
		// gsl_sf_bessel_K_scaled_steed_temme_CF2 x > 2
		result[Dcol * sizeResultPadRow + Drow] = logthisx;
		if(logthisx > 2.0) { // gsl_sf_bessel_K_scaled_steed_temme_CF2
			K_nu = 0.0;
			bi = 2.0*(1.0 + thisx);//x);"
			di = 1.0/bi;

			delhi = di;
			hi = delhi;// was di, now di/exp(maternBit)

			qi   = 0.0;
			qip1 = 1.0;

			ai = -(0.25 - muSq);
			a1 = ai;
			ci = -ai;
			Qi = -ai;

			s = 1.0 + Qi*delhi;

			for(k=2; k<=maxIter; k++) {
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
				if(fabs(dels/s) < epsilon) break;
				} // k loop

			hi *= -a1;
			K_nu = exp(-thisx) * exp(maternBit) * sqrt(M_PI/(2.0*thisx)) / s;//"  sqrt(pi)/2 sqrt(2/x)/s =

			K_nup1 = K_nu * (mu + thisx + 0.5 - hi)/thisx;
			Kp_nu  = - K_nup1 + mu/thisx * K_nu;

			} else { // if short distance gsl_sf_bessel_K_scaled_temme
			twoLnHalfX = 2*ln_half_x;
			sigma   = - mu * ln_half_x;
			half_x_nu = exp(-sigma);
			sinhrat = sinh(sigma)/sigma;

			k = sinrat * (cosh(sigma)*g1 - sinhrat*ln_half_x*g2);
			pk = 0.5/half_x_nu * g_1pnu;
			qk = 0.5*half_x_nu * g_1mnu;
			hk = pk;
			sum0 = fk*expMaternBit;
			sum1 = hk*expMaternBit;
			k=0;
			logck = maternBit;
			del0 = fabs(sum0)+100;
			ck, del1;

			while( (k < maxIter) && ( fabs(del0) > (epsilon * fabs(sum0)) ) ) {
				k++;
				logck += twoLnHalfX - log((double)k);
				ck = exp(logck);
				fk  = (k*fk + pk + qk)/(k*k-muSq);

				del0 = ck * fk;
				sum0 += del0;

				pk /= (k - (mu));
				qk /= (k + (mu));

				hk  = -k*fk + pk;
				del1 = ck * hk; 	
				sum1 += del1;
				}//while loop
			//		result[Dcol * sizeResultPadRow + Drow] = -k;
			K_nu   = sum0;
			K_nup1 = sum1 * exp( - ln_half_x);
			} // short distance

		for(k=0; k<nuround; k++) {
			K_num1 = K_nu;
			K_nu   = K_nup1;
			K_nup1 = exp(log(mup1+k) - ln_half_x) * K_nu + K_num1;
		}
		result[Dcol + Drow * sizeResultPadRow] = K_nu; // lower triangle
		result[Drow + Dcol * sizeResultPadRow] = K_nu; // upper triangle
		}// not diag

};//function
