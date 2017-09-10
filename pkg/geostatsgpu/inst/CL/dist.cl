__kernel void distCL(
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
		__global double *result)
{
	// Get the index of the elements to be processed
	const int Drow = get_global_id(0);
	const int Dcol = get_global_id(1);

	if(Drow < size && Dcol < Drow){
		result[Dcol + Drow * sizeResultPadCol] = Dcol;
		result[Dcol * sizeResultPadRow + Drow] = Drow;
	}//not diag

};//function
