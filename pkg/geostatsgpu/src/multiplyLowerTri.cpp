#include "geostatsgpu.hpp"



// C = A B, A lower triangular

template <typename T> 
std::string multiplyLowerBatchString(
  const int sameB,
  const int Nrow, 
  const int Ncol,
  const int Nmatrix, 
  const int NpadC, 
  const int NpadA, 
  const int NpadB,
  const int NpadBetweenMatricesC,
  const int NpadBetweenMatricesA,
  const int NpadBetweenMatricesB,
  const int NlocalCacheA,  // greater than Nlocal(0), smaller than Nrow
  const int NlocalCacheB) {

  std::string typeString = openclTypeString<T>();
  std::string result = "";

  if(typeString == "double") {
  	result += "\n#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n\n";
  }
  result = result + 
"\n#define Nrow " + std::to_string(Nrow) + "\n"    
"#define Ncol " + std::to_string(Ncol) + "\n"    
"#define Nmatrix " + std::to_string(Nmatrix) + "\n"    
"#define NpadC " + std::to_string(NpadC) + "\n"    
"#define NpadA " + std::to_string(NpadA) + "\n"    
"#define NpadB " + std::to_string(NpadB) + "\n"    
"#define NpadBetweenMatricesC " + std::to_string(NpadBetweenMatricesC) + "\n"    
"#define NpadBetweenMatricesA " + std::to_string(NpadBetweenMatricesA) + "\n"    
"#define NpadBetweenMatricesB " + std::to_string(NpadBetweenMatricesB) + "\n"    
"#define NlocalCache " + std::to_string(NlocalCacheB) + "\n"    
"#define NlocalCacheA " + std::to_string(NlocalCacheA) + "\n\n"    

// global work items are rows, columns, matrices
// local work items anything, anything, 1
"__kernel void multiplyLowerBatch(\n"
  "	__global " + typeString+ " *C,\n"
  "	__global "+ typeString+ " *A,\n"
  "	__global "+ typeString+ " *B) {\n\n" +

    typeString + " Dout, *AHere, *BHere, *CHere;\n"
  "	local "+ typeString+ " Acache[NlocalCacheA];\n" 
"	local "+ typeString+ " Bcache[NlocalCache];\n"
  
"int Dmatrix, Drow, Dcol, Dinner, DrowNpadC, DrowNpadA;\n"
  
"const int doCacheA = (get_local_id(1) == 0);\n"
"const int DinnerStop = min(Nrow, NlocalCache);\n";

  "BHere = B;\n"

    "for(Dmatrix = get_global_id(2); Dmatrix < Nmatrix; Dcol += get_global_size(2)){\n"
  "  AHere = &A[Dmatrix*NpadBetweenMatricesA];\n";

  if(!sameB){
    result +=  "  BHere = &B[Dmatrix*NpadBetweenMatricesB];\n";
  } 

result += "  CHere = &C[Dmatrix*NpadBetweenMatricesC];\n"

// cache first rows of B 
"  for(Drow = get_local_id(0); Drow < NlocalCache; Drow += get_local_size(0)){\n"
"    DrowNpadC = Drow*Ncol;\n"
"    DrowNpadA = Drow*NpadB;\n"
"    for(Dcol = get_local_id(1); Dcol < Ncol; Dcol += get_local_size(1))){\n"
//     "Bcache[Dcol + Drow*Ncol] = B[Dcol + Drow * NpadB];\n"
"      Bcache[Dcol + DrowNpadC] = BHere[Dcol +DrowNpadA];\n"
"    }\n"
"  }\n"
"barrier(CLK_LOCAL_MEM_FENCE);\n"

// looped through rows which are all cached
"for(Drow = get_global_id(0); Drow < DinnerStop; Drow+=get_global_size(0)){\n"
  "DrowNpadA= Drow * NpadA;\n"
  "DrowNpadC= Drow * NpadC;\n"
  "for(Dcol = get_global_id(1); Dcol < Ncol; Dcol += get_global_size(1)){\n"
    "Dout = 0.0;\n"
    "for(Dinner = 0; Dinner <= Drow; Dinner++){\n"
    "  if(doCacheA) {\n"
    "    Acache[get_local_id(0)] = AHere[Dinner + DrowNpadA];\n"
    "  }\n"
	  "  barrier(CLK_LOCAL_MEM_FENCE);\n"
	  "Dout += Acache[get_local_id(0)] * Bcache[Dcol + Dinner * Ncol];\n"
	"}\n" // Dinner
  	"CHere[Dcol + DrowNpadC] = Dout;\n"
  "}\n" // Dcol
"}\n" //Drow
// rows which are not all cached
"for( ; Drow < Nrow; Drow+=get_global_size(0)){\n"
  "DrowNpadA= Drow * NpadA;\n"
  "DrowNpadC= Drow * NpadC;\n"
  "for(Dcol = get_global_id(1); Dcol < Ncol; Dcol += get_global_size(1)){\n"
    "Dout = 0.0;\n"
    // cached rows
    "for(Dinner = 0; Dinner < NlocalCache; Dinner++){\n"
	  "if(doCacheA) {\n"
	    "Acache[get_local_id(0)] = AHere[Dinner + DrowNpadA];\n"
	  "}\n"
	  "barrier(CLK_LOCAL_MEM_FENCE);\n"
	  "Dout += Acache[get_local_id(0)] * Bcache[Dcol + Dinner * Ncol];\n"
	"}\n" // Dinner
	// un-cached rows
    "for( ; Dinner <= Drow; Dinner++){\n"
	  "if(doCacheA) {\n"
	    "Acache[get_local_id(0)] = AHere[Dinner + DrowNpadA];\n"
	  "}\n"
	  "barrier(CLK_LOCAL_MEM_FENCE);\n"
	  "Dout += Acache[get_local_id(0)] * BHere[Dcol + Dinner * NpadB];\n"
	"}\n" // Dinner
	"CHere[Dcol + DrowNpadC] = Dout;\n"
  "}\n" // Dcol
 "}\n" //Drow
 
 "}\n" // Dmatrix
 "}";

  return(result);
}





template <typename T> 
void multiplyLowerBatch(
	viennacl::matrix<T> &C,
	viennacl::matrix<T> &A,
	viennacl::matrix<T> &B,
	std::vector<int> Nglobal, 
	std::vector<int> Nlocal, 
	const int NlocalCache, 
	const int ctx_id) {

  const int Nmatrix = A.size1()/A.size2();
  const int sameB = (A.size2() == B.size1());
	// the context
	viennacl::ocl::context ctx(viennacl::ocl::get_context(ctx_id));

	cl_device_type type_check = ctx.current_device().type();

	viennacl::ocl::program & my_prog = ctx.add_program(
	  multiplyLowerBatchString<T>(  
	    sameB,
  	  A.size2(), 
	    B.size2(),
	    Nmatrix, 
  	  C.internal_size2(), 
	    A.internal_size2(), 
	    B.internal_size2(),
  	  C.internal_size2()*Nmatrix,//NpadBetweenMatricesC,
	    A.internal_size2()*Nmatrix,//NpadBetweenMatricesA,
	    B.internal_size2()*Nmatrix,//NpadBetweenMatricesB,
  	  Nlocal[0], 
	    NlocalCache),
	    "my_kernel");
	
	viennacl::ocl::kernel & multiplyLowerKernel = my_prog.get_kernel("multiplyLowerBatch");

  multiplyLowerKernel.global_work_size(0, Nglobal[0]);
  multiplyLowerKernel.global_work_size(1, Nglobal[1]);
  multiplyLowerKernel.global_work_size(2, Nglobal[2]);

  multiplyLowerKernel.local_work_size(0, Nlocal[0]);
  multiplyLowerKernel.local_work_size(1, Nlocal[1]);
  multiplyLowerKernel.local_work_size(2, Nlocal[2]);
  
		// diagonals and diagTimesRowOfA
		viennacl::ocl::enqueue(multiplyLowerKernel(
			C, A, B));
		
}

template <typename T> 
SEXP multiplyLowerBatchTyped(
	Rcpp::S4 CR,
	Rcpp::S4 AR,
	Rcpp::S4 BR,
	Rcpp::IntegerVector NglobalR,
	Rcpp::IntegerVector NlocalR, 
	int NlocalCache) {

  std::vector<int> Nglobal = Rcpp::as<std::vector<int> >(NglobalR);
  std::vector<int> Nlocal = Rcpp::as<std::vector<int> >(NlocalR);

	const int ctx_id = INTEGER(CR.slot(".context_index"))[0]-1;
	const bool BisVCL=1;


	
	std::shared_ptr<viennacl::matrix<T> > 
		AG = getVCLptr<T>(AR.slot("address"), BisVCL, ctx_id);
	std::shared_ptr<viennacl::matrix<T> > 
		BG = getVCLptr<T>(BR.slot("address"), BisVCL, ctx_id);
	std::shared_ptr<viennacl::matrix<T> > 
		CG = getVCLptr<T>(CR.slot("address"), BisVCL, ctx_id);

	multiplyLowerBatch<T>(*CG, *AG, *BG, Nglobal, Nlocal, NlocalCache, ctx_id);	

	return Rcpp::wrap(0L);
}


//' Multiply lower triangular matrices
//' 
//' Multiplies a lower triangular matrix by a rectangular matrix
//'
//' @param C output matrices, stacked row-wise
//' @param A lower triangular matrices
//' @param B rectangular matrix or matrices
//' @param Nglobal vector of number of global work items
//' @param Nlocal vector of number of local work items
//' @param NlocalCache elements in local cache
//' @export
// [[Rcpp::export]]
SEXP multiplyLowerBatchBackend(
	Rcpp::S4 C,
	Rcpp::S4 A,
	Rcpp::S4 B,
	Rcpp::IntegerVector Nglobal,
	Rcpp::IntegerVector Nlocal, 
	int NlocalCache) {

	SEXP result;
  
  Rcpp::traits::input_parameter< std::string >::type classVarR(RCPP_GET_CLASS(C));
  std::string precision_type = (std::string) classVarR;
  

  if(precision_type == "fvclMatrix") {
    result = multiplyLowerBatchTyped<float>(C, A, B, Nglobal, Nlocal, NlocalCache);
  } else if (precision_type == "dvclMatrix") {
    result = multiplyLowerBatchTyped<double>(C, A, B, Nglobal, Nlocal, NlocalCache);
	} else {
		result = Rcpp::wrap(1L);
	}
	return(result);

}
