/*
 * xScale = x * sqrt(8 * shape) / range
 * xScale = x * 2 * sqrt(2 * shape) / range
 * log(xscale) = log(x) + log(2) + 0.5*log(2) + 
 *   0.5 * log(shape) - log(range)
 * log(xscale) = logxscale + (1/2) * log(x^2)
 * logxscale = 1.5*log(2) + 0.5 * log(shape) - log(range)
 * matern(x; sigma, range, shape) = 
 *   (sigma^2 / Gamma(shape)) * 2^(shape-1) * xScale^ shape *  
 *   besselK(xScale, shape)
 */

__kernel void maternCL(
		const int size,
		const int sizeCoordsPadCol,
		const int sizeResultPadRow,
		const int sizeResultPadCol,
		const double nu,
		const int nuround,
		const double mu,
		const double costheta,
		const double sintheta,
		const double anisoRatioSq,
		const double varscale,
		const double logxscale,
		const double diagVar,
		const double sinrat,
		const double g_1pnu,
		const double g_1mnu,
		const double g1,
		const double g2,
		const double epsilon,
		__global const double *coords,
		__global double *result) {
	// Get the index of the elements to be processed
	const int Drow = get_global_id(0);
	const int Dcol = get_global_id(1);
	const int local0 = get_local_id(0);
	const int local1 = get_local_id(1);
	if(Drow < size && Dcol == Drow) {
		result[Dcol + Drow * sizeResultPadRow] = diagVar; // lower triangle

	}
	if(Drow < size && Dcol < Drow){

		double dist[2], distRotate[2];
		dist[0] = coords[Drow*sizeCoordsPadCol] - coords[Dcol*sizeCoordsPadCol];
		dist[1] = coords[Dcol*sizeCoordsPadCol +1] - coords[Drow*sizeCoordsPadCol +1];
		distRotate[0] = costheta *dist[0] - sintheta * dist[1];
		distRotate[1] = sintheta *dist[0] + costheta * dist[1];
		const double distSq = distRotate[0]*distRotate[0] + distRotate[1]*distRotate[1]/anisoRatioSq;


		const double logthisx = log(distSq)/2 + logxscale;
		double K_mu, K_mup1, logck;
		double K_nu=-1.1, K_nup1, K_num1, Kp_nu;
		int k;
		double twoLnHalfX, sigma, half_x_nu, sinhrat;

		const double ln_half_x = logthisx - M_LN2;
		const double thisx = exp(logthisx);
		const double maternBit = varscale + nu * logthisx;
		const double expMaternBit = exp(maternBit);

		const int maxIter = 15000;
		const double muSq = mu * mu, mup1 = mu+1;

		// gsl_sf_bessel_K_scaled_temme x < 2
		// gsl_sf_bessel_K_scaled_steed_temme_CF2 x > 2
//		result[Dcol * sizeResultPadRow + Drow] = logthisx;
		if(logthisx > 2.0) { // gsl_sf_bessel_K_scaled_steed_temme_CF2
			K_nu = 0.0;
			double bi = 2.0*(1.0 + thisx);//x);"
//			double logbi = M_LN2 + log1p(thisx);
//			double bi = exp(logbi);
//			double di = exp(-logbi);//1.0/bi;
			double di = 1.0/bi;
			//		  "double delhi = di;" // divide by exp(maternBit)
//			double delhi = exp(-logbi-maternBit);
			double delhi = di;
			double hi    = delhi;// was di, now di/exp(maternBit)

			double qi   = 0.0;
			double qip1 = 1.0;

			double ai = -(0.25 - muSq);
			double a1 = ai;
			double ci = -ai;
			double Qi = -ai;

//			double s = exp(-maternBit) + Qi*delhi;
			double s = 1.0 + Qi*delhi;
			double dels;
			double tmp;
			for(k=2; k<=maxIter; k++) {
//				dels = 0.0;
//				tmp = 0.0;
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
			//		result[Dcol * sizeResultPadRow + Drow] = -k;

			hi *= -a1;

//			K_nu= exp(-0.5*ln_half_x)/(M_2_SQRTPI * s);
			K_nu   = exp(-thisx) * exp(maternBit) * sqrt(M_PI/(2.0*thisx)) / s;//"  sqrt(pi)/2 sqrt(2/x)/s =

//			K_nup1 = K_nu * (mu + thisx + 0.5 - hi*exp(maternBit))/thisx;
			K_nup1 = K_nu * (mu + thisx + 0.5 - hi)/thisx;
			Kp_nu  = - K_nup1 + mu/thisx * K_nu;

		} else { // if short distance gsl_sf_bessel_K_scaled_temme
			twoLnHalfX = 2*ln_half_x;
			sigma   = - mu * ln_half_x;
			half_x_nu = exp(-sigma);
			sinhrat = sinh(sigma)/sigma;

			double fk = sinrat * (cosh(sigma)*g1 - sinhrat*ln_half_x*g2);
			double pk = 0.5/half_x_nu * g_1pnu;
			double qk = 0.5*half_x_nu * g_1mnu;
			double hk = pk;
			double sum0 = fk*expMaternBit;
			double sum1 = hk*expMaternBit;
			k=0;
			double logck = maternBit;
			double del0 = fabs(sum0)+100;
			double ck, del1;

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
//		result[Dcol * sizeResultPadRow + Drow] = logthisx; // upper triangle

	}// not diag

};//function
