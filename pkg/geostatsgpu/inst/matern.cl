static const char * kernel_matern =
		__kernel void matern(
				const int size,
				const int internalSizeCoords,
				const int internalSizeResult,
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
				__global double *result)
{
	// Get the index of the elements to be processed
	const int globalRow = get_global_id(0);
	if(globalRow < size){
		const int rowHereResult = internalSizeResult*globalRow;
		const int rowHereCoords1 = internalSizeCoords*globalRow;
		int Dcol, rowHereCoords2;
		const int maxIter = 15000;
		const double muSq = mu * mu, mup1 = mu+1;
		double dist[2], distRotate[2], distSq;
		double del0, del1, sum0, sum1,fk, pk, qk, hk, ck;
		double K_mu, K_mup1, logck;
		double K_nu, K_nup1, K_num1, Kp_nu;
		int k;
		double logthisx, ln_half_x, twoLnHalfX, maternBit, sigma, half_x_nu, sinhrat;

		result[globalRow+rowHereResult] = diagVar;
		for(Dcol=0; Dcol < globalRow; ++Dcol){
			rowHereCoords2 = internalSizeCoords*Dcol;
			dist[0] = coords[rowHereCoords1] - coords[rowHereCoords2];
			dist[1] = coords[rowHereCoords1 + 1] - coords[rowHereCoords2 + 1];
			distRotate[0] = costheta *dist[0] - sintheta * dist[1];
			distRotate[1] = sintheta *dist[0] + costheta * dist[1];
			distSq = distRotate[0]*distRotate[0] + distRotate[1]*distRotate[1]/anisoRatioSq;

			logthisx = 0.5 * log(distSq) + logxscale;
			ln_half_x = logthisx - M_LN2;
			twoLnHalfX = 2*ln_half_x;
			maternBit = varscale + nu * logthisx;
			sigma   = - mu * ln_half_x;
			half_x_nu = exp(-sigma);
			sinhrat = sinh(sigma)/sigma;

			fk = sinrat * (cosh(sigma)*g1 - sinhrat*ln_half_x*g2);
			pk = 0.5/half_x_nu * g_1pnu;
			qk = 0.5*half_x_nu * g_1mnu;
			hk = pk;
			sum0 = fk*exp(maternBit);
			sum1 = hk*exp(maternBit);
			k=0;
			logck = maternBit;
			del0 = fabs(sum0)+100;

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

			K_nu   = sum0;
			K_nup1 = sum1 * exp( - ln_half_x);

			for(k=0; k<nuround; k++) {
				K_num1 = K_nu;
				K_nu   = K_nup1;
				K_nup1 = exp(log(mup1+k) - ln_half_x) * K_nu + K_num1;
			}
			result[Dcol+rowHereResult] = K_nu;
		} // loop col
	} // if size
}//function
