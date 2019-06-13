#include "geostatsgpu.hpp"
#define DEBUG

template <typename T> 
int sizeOfReal() {
  return(-1);}

template <> int sizeOfReal<double>(){
  return(sizeof(cl_double));}

template <> int sizeOfReal<float>(){
  return(sizeof(cl_float));}



template <typename T> 
std::string openclTypeString() {
  return("undefined");}

template <> std::string openclTypeString<double>(){
  std::string result = "double";
  return(result);
}

template <> std::string openclTypeString<float>(){
  std::string result = "float";
  return(result);
}

template <> std::string openclTypeString<int>(){
  std::string result = "int";
  return(result);
}



template <typename T> 
std::string maternBatchKernelString(int maxIter) {

  std::string typeString = openclTypeString<T>();
  std::string maxIterString = std::to_string(maxIter);
  std::string result = "";

  if(typeString == "double") {
  	result += "\n#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n\n"
  "\n#define logSqrtHalfPi 0.2257913526447273278031\n"
  "\n#define M_PI_T M_PI\n"
  "\n#define M_LN2_T M_LN2\n";
  "\n#define epsilon " + std::to_string(GSL_DBL_EPSILON/1000) + "\n";
  } else {
  	result += "\n#define logSqrtHalfPi 0.22579135\n"
	"\n#define M_LN2_T M_LN2_F\n"
    "\n#define M_PI_T M_PI_F\n"
	"\n#define epsilon " + std::to_string(GSL_FLT_EPSILON /1000) + "\n";
  }

  result += "\n#define maxIter " + maxIterString + "\n"
	"#define assignDiag\n"
	"#define assignLower\n"
	"#define assignLower\n";




  result = result + 
"\n void maternShort(\n" +
typeString + " ln_half_x, " + 
typeString + " maternBit, " + typeString + " expMaternBit, " + 
typeString + " mu, " + typeString + " muSq, " + 
typeString + " sinrat, " +
typeString + " g1, " + typeString + " g2," +
typeString + " g_1pnu, " + typeString + " g_1mnu," +
typeString + " *K_nu, " + typeString + " *K_nup1\n" 
"  ){\n" +
typeString + " twoLnHalfX, sigma, sinhrat, fk, pk, qk, hk, half_x_nu;\n" +
typeString + " sum0, sum1, ck, logck, del0, del1;\n"
"int k;\n"
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
	"			*K_nu   = sum0;\n"
	"			*K_nup1 = sum1 * exp( - ln_half_x);\n"
"}\n\n"


"\n void maternLong(\n" +
typeString + " thisx, " + typeString + " expMaternBit, " + 
typeString + " nu, " + typeString + " mu, " + typeString + " muSq, " + 
typeString + " *K_nu, " + typeString + " *K_nup1\n"
"){\n" +
typeString + " bi, di, delhi, hi, qi, qip1, ai, ci, Qi, s, tmp;\n" +
typeString + " dels, a1;\n"
"int k;\n"
	"		*K_nu = 0.0;\n"
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
	"		*K_nu = exp(-thisx) * expMaternBit * sqrt(M_PI_T/(2.0*thisx)) / s;\n"//  sqrt(pi)/2 sqrt(2/x)/s =

	"		*K_nup1 = *K_nu * (mu + thisx + 0.5 - hi)/thisx;\n"
// not needed?	"		Kp_nu  = - K_nup1 + mu/thisx * K_nu;\n"
"}\n\n"

// parameters
// 0 range, 1 shape, 2 variance, 3 nugget, 
// 4 anisoRatio 5 anisoAngleRadians  6 anisoAngleDegrees
// 7 sintheta, 8 costheta 9 anisoRatioSq
// 10 varscale
// 11 logxscale
// 12 sinrat 13 mu 14 muSq 15 mup1 16 nuround
// 17 g1 18 g2 19 g1pnu 20 g1mnu
// 21 variance + nugget
"\n#define NlocalParams 22\n"
"\n__kernel void maternBatch(\n"
  "__global " + typeString + " *result,"
  "__global " + typeString + " *coords,\n"
  "const int Ncell, const int N, const int Npad," // padding for inner matrices, usually N*N
  "const int Nmatrix, const int NpadResult, const int NpadCoords,\n"
  "__global " + typeString + " *params, const int NpadParams,\n"  
  "__local " + typeString + " *localParams,\n"
  "__local " + typeString + " *localCoords"
  ") {\n"

"__local int Drow, Dcol;\n"

"int Dmatrix, Dcell, nuround, DlocalParam, k;\n" +
typeString + " distRotate[2], distSq;\n" +
typeString + " logthisx, ln_half_x, thisx, maternBit, expMaternBit;\n" +
typeString + " K_num1, K_nu, K_nup1;\n\n"

// dimension 0 is cell, dimension 1 is matrix
"int isFirstLocal = (get_local_id(0)==0 & get_local_id(1)==0);\n\n"
"int isFirstLocal1 = (get_local_id(1)==0);\n\n"
"int distanceY = get_local_id(0) + get_local_size(0);\n\n\n"

// copy parameters to local storage
"if(get_local_id(0)==0){\n"
"for(Dmatrix = get_local_id(1); Dmatrix < Nmatrix; Dmatrix += get_local_size(1)) {\n"
	"DlocalParam = NlocalParams*Dmatrix;\n"
	"k = NpadParams*Dmatrix;\n"

	"for(Dcell = 0; Dcell < N; ++Dcell){\n"
		"localParams[DlocalParam + Dcell] = params[k + Dcell];\n"
	"}\n" // Dcell
"}\n" // Dmatrix
"}\n" // if local0


// copy coordinates to local storage
"if(get_local_id(0)==1){\n"
"for(Dmatrix = get_local_id(1); Dmatrix < Nmatrix; Dmatrix += get_local_size(1)) {\n"
	"DlocalParam = NpadCoords*Dmatrix ;\n"
	"k = 2*Dmatrix;\n"
	"for(Dcell = 0; Dcell < N; ++Dcell){\n"
		"localCoords[k+ Dcell] = coords[DlocalParam + Dcell];\n"
	"}\n" // Dcell
"}\n" // Dmatrix
"}\n" // if local0
"barrier(CLK_LOCAL_MEM_FENCE);\n"
"for(Dcell = get_global_id(0); Dcell < Ncell; Dcell += get_global_size(0)) {\n"


"if(isFirstLocal){\n"
	"	Drow = ceil(0.5 + sqrt(0.25 + 2.0*(Dcell+1) ) ) - 1;\n"
	"	Dcol = Dcell - round(Drow * (Drow - 1.0) / 2.0);\n"
"}\n\n" // if local0
"barrier(CLK_LOCAL_MEM_FENCE);\n"

"if(isFirstLocal1){\n"

	"	localCoords[get_local_id(0)] = coords[Drow*NpadCoords] - coords[Dcol*NpadCoords];\n"
	"	localCoords[distanceY] = coords[Dcol*NpadCoords +1] - coords[Drow*NpadCoords +1];\n"
"}\n\n" 
"barrier(CLK_LOCAL_MEM_FENCE);\n"
#ifdef UNDEF

"for(Dmatrix = get_global_id(1); Dmatrix < Nmatrix; Dcell += get_global_size(1)) {\n"
	"   DlocalParam = NlocalParams*Dmatrix;\n"
	"   nuround = (int) (localParams[DlocalParam+16]);\n"

	"	distRotate[0] = localParams[DlocalParam+8] * localCoords[get_local_id(0)] +"
	"         localParams[DlocalParam+7] * localCoords[distanceY];\n"
	"	distRotate[1] = localParams[DlocalParam+7] * localCoords[get_local_id(0)] -" 
	"         localParams[DlocalParam+8] * localCoords[distanceY];\n"

	"	distSq = distRotate[0]*distRotate[0] +" 
	"         distRotate[1]*distRotate[1]/localParams[DlocalParam + 9];\n"

	"	logthisx = log(distSq)/2 + localParams[DlocalParam + 11];\n"
	"	ln_half_x = logthisx - M_LN2_T;\n"
	"	thisx = exp(logthisx);\n"
	"	maternBit = localParams[DlocalParam + 10] + localParams[DlocalParam + 1] * logthisx;\n"
	"	expMaternBit = exp(maternBit);\n"

	"	if(logthisx > 2.0) {\n" // gsl_sf_bessel_K_scaled_steed_temme_CF2

//	"   maternLong(thisx, expMaternBit, nu[Dmatrix], mu[Dmatrix], muSq[Dmatrix], &K_nu, &K_nup1);\n"
	"   maternLong(thisx, expMaternBit, localParams[DlocalParam + 1],"
	"       localParams[DlocalParam + 13], localParams[DlocalParam + 14], &K_nu, &K_nup1);\n"

	"	} else { \n"// if short distance gsl_sf_bessel_K_scaled_temme
			"maternShort(ln_half_x, maternBit, expMaternBit, "
			//" mu[Dmatrix], muSq[Dmatrix]," 
			" localParams[DlocalParam + 13], localParams[DlocalParam + 14],"
			" localParams[DlocalParam + 12]," // sinrat
//			" g1[Dmatrix], g2[Dmatrix],"
			" localParams[DlocalParam + 17], localParams[DlocalParam + 18],"
//			" g_1pnu[Dmatrix], g_1mnu[Dmatrix]," 
			" localParams[DlocalParam + 19], localParams[DlocalParam + 20],"
			" &K_nu, &K_nup1);"
	"   }\n"

	"		for(k=0; k<nuround; k++) {\n"
	"			K_num1 = K_nu;\n"
	"			K_nu   = K_nup1;\n"
	"			K_nup1 = exp(log(mup1[Dmatrix]+k) - ln_half_x) * K_nu + K_num1;\n"
	"		}\n"
"#ifdef assignLower\n"
"	  result[Dcol + Drow * sizeResultPadRow] = K_nu;\n" // lower triangle
"# endif\n"
"#ifdef assignUpper\n;"
"	result[Drow + Dcol * sizeResultPadRow] = K_nu;\n"//K_nu;\n" // upper triangle
"# endif\n"



	"}\n" // Dmatrix
#endif

"}\n" // Dcell
//"\n#ifdef assignDiag\n;"


"barrier(CLK_LOCAL_MEM_FENCE);\n"

"for(Dmatrix = get_global_id(1); Dmatrix < Nmatrix; Dmatrix += get_global_size(1)) {\n"
	"DlocalParam = Dmatrix * NpadResult;\n"
	"maternBit = localParams[NlocalParams*Dmatrix+21];\n"
	///FIX THIS!!!
//	"maternBit = params[NpadParams*Dmatrix+21];\n"
	"for(Dcell = get_global_id(0); Dcell < N; Dcell += get_global_size(0)) {\n"
	"	result[Dcell + Dcell * Npad + DlocalParam] = maternBit;\n"
	"}\n" // Dmatrix

"}\n" // Dcell
//"\n#endif\n"


"}\n"; // function
return result;
}

template<typename T> 
void fill22params(
	viennacl::matrix<T> &param
	// Nmat rows, 22 columns
// 0 range, 1 shape, 2 variance, 3 nugget, 
// 4 anisoRatio 5 anisoAngleRadians  6 anisoAngleDegrees
// 7 sintheta, 8 costheta 9 anisoRatioSq
// 10 varscale
// 11 logxscale
// 12 sinrat 13 mu 14 muSq 15 mup1 16 nuround
// 17 g1 18 g2 19 g1pnu 20 g1mnu
// 21 variance + nugget

) {

	int D, nuround;
	const int Nmat = param.size1();
	T range, shape, theta, anisoRatio, variance;
	T onePointFiveM_LN2 = 1.5 * M_LN2;
	double muSq, pi_nu, mu;
	double g_1pnu, g_1mnu, g1, g2;
	T epsHere = maternClEpsilon<T>();


	for(D=0; D<Nmat;++D) {

		range = param(D, 0);
		shape = param(D, 1);
		variance = param(D, 2);
		anisoRatio = param(D, 4);
		theta = param(D, 5);

		nuround = (int) (shape+0.5);
		mu = shape - nuround;
		pi_nu = M_PI * mu,

		Rtemme_gamma(&mu, &g_1pnu, &g_1mnu, &g1, &g2);


		param(D, 7) = (T) cos(theta);
		param(D, 8) = (T) sin(theta);
		param(D, 9) = (T) anisoRatio*anisoRatio;
		// varscale
		param(D, 10) = (T) (log(variance) - Rf_lgammafn(shape) - (shape-1)*M_LN2);
		// logxscale
		param(D, 11) = (T) (onePointFiveM_LN2 + 0.5 * log(shape) - log(range));
		// sinrat
		param(D, 12) = (T) (fabs(pi_nu) < epsHere ? 1.0 : pi_nu/sin(pi_nu));
		// mu 
		param(D, 13)= (T) mu;		
		// muSq
		param(D, 14)= (T) mu*mu;		
		// mup1 
		param(D, 15)= (T) mu + 1.0;		
		// nuround
		param(D, 16)= (T) nuround;		
		// g1 
		param(D, 17)= (T) g1;		
		// g2 
		param(D, 18)= (T) g2;		
		// g1pnu 
		param(D, 19)= (T) g_1pnu;		
		// g1mnu
		param(D, 20)= (T) g_1mnu;		
		// variance + nugget
		param(D, 21)= (T) variance + param(D, 3);		
	}

}

template<typename T> void maternBatchVcl(
	viennacl::matrix<T> &vclVar, // Nmat columns N^2 rows
	viennacl::matrix<T> &vclCoords, // 2 columns
	viennacl::matrix<T> &param, // Nmat rows, 22 columns
	std::vector<int> numWorkItems,
	std::vector<int> numLocalItems,	
	const int ctx_id)
{

	const int 
	N = vclCoords.size1(),
	Nmat = vclVar.size2(),
	iSizeCoords2=vclCoords.internal_size2(),
	iSizeVar2=vclVar.internal_size2(),
	iSizeParam2=param.internal_size2();

	const int Ncell = N * (N - 1)/2, maxIter = 1500;
	const int Npad = N; // padding for matrix rows, can't be deduced from VCL object

	fill22params(param);

	// the context
	viennacl::ocl::context ctx(viennacl::ocl::get_context(ctx_id));

	cl_device_type type_check = ctx.current_device().type();

	std::string maternClString = maternBatchKernelString<T>(maxIter);



	viennacl::ocl::program & my_prog = ctx.add_program(maternClString, "my_kernel");
	// get compiled kernel function
	viennacl::ocl::kernel & maternKernel = my_prog.get_kernel("maternBatch");

	// dimension 0 is cell, dimension 1 is matrix
	maternKernel.global_work_size(0, (cl_uint) (numWorkItems[0] ) );//numWorkItems[0]);
	maternKernel.global_work_size(1, (cl_uint) (numWorkItems[1] ) );//numWorkItems[0]);

	maternKernel.local_work_size(0, (cl_uint) (numLocalItems[0]));
	maternKernel.local_work_size(1, (cl_uint) (numLocalItems[1]));

	// local memory
	viennacl::ocl::local_mem localParams(
		param.size1() * param.size2() * sizeOfReal<T>());
	viennacl::ocl::local_mem localCoords(
		vclCoords.size1() * vclCoords.size2() * sizeOfReal<T>());

	Rcpp::Rcout << maternClString;

#ifdef DEBUG

Rcpp::Rcout << "Ncell " << Ncell << " N" << N <<" Npad" << Npad << " Nmatrix" << Nmat << " NpadResult" << iSizeVar2  << " NpadCoords" << iSizeCoords2  << " NpadParams" << iSizeParam2;

#endif
	viennacl::ocl::enqueue(maternKernel(
		vclVar, vclCoords,
// Ncell,  Npad  Nmatrix  NpadResult  NpadCoords params NpadParams
  Ncell, N, Npad, Nmat, iSizeVar2, iSizeCoords2, 
		param, iSizeParam2, 
		localParams, localCoords 
	));


}


template<typename T> void maternBatchTemplated(
	Rcpp::S4 varR,
	Rcpp::S4 coordsR,
	Rcpp::S4 paramR, //22 columns 
	Rcpp::IntegerVector numWorkItems,
	Rcpp::IntegerVector numLocalItems
) {

	std::vector<int> numWorkItemsStd = Rcpp::as<std::vector<int> >(numWorkItems);
	std::vector<int> numLocalItemsStd = Rcpp::as<std::vector<int> >(numLocalItems);
	
	// data
	const bool BisVCL=1;
	const int ctx_id = INTEGER(varR.slot(".context_index"))[0]-1;
	std::shared_ptr<viennacl::matrix<T> > vclVar = getVCLptr<T>(varR.slot("address"), BisVCL, ctx_id);
	std::shared_ptr<viennacl::matrix<T> > param = getVCLptr<T>(paramR.slot("address"), BisVCL, ctx_id);
	std::shared_ptr<viennacl::matrix<T> > vclCoords = getVCLptr<T>(coordsR.slot("address"), BisVCL, ctx_id);

	maternBatchVcl<T>(
		*vclVar, *vclCoords,
		*param,
        numWorkItemsStd, 
        numLocalItemsStd,
        ctx_id);

}

//[[Rcpp::export]]
void maternBatchBackend(
	Rcpp::S4 varR,
	Rcpp::S4 coordsR,
	Rcpp::S4 paramR, //22 columns 
	Rcpp::IntegerVector numWorkItems,
	Rcpp::IntegerVector numLocalItems) {


    Rcpp::traits::input_parameter< std::string >::type classVarR(RCPP_GET_CLASS(varR));
    std::string precision_type = (std::string) classVarR;


  if(precision_type == "fvclMatrix") {
	maternBatchTemplated<float>(
		varR, coordsR, paramR, 
		numWorkItems, numLocalItems);
  } else if (precision_type == "dvclMatrix") {
	maternBatchTemplated<double>(
		varR, coordsR, paramR, 
		numWorkItems, numLocalItems);
  } else {
  	Rcpp::warning("precision_type must be float or double");
 }

}
