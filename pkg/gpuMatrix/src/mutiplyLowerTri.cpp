

#include "gpuMatrix.hpp"

// C = A B, A lower triangular


template <typename T> 
std::string multiplyLowerTypeString() {
	
  std::string typeString = openclTypeString<T>();
  std::string result = "";

  if(typeString == "double") {
  	result += "\n#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n\n";
  }
  result = result + 
"__kernel void multiplyLower(\n"
  "	__global " + typeString+ " *C,\n"
  "	__global "+ typeString+ " *A,\n"
  "	__global "+ typeString+ " *B,\n"
 "	__local "+ typeString+ " *Acache,\n" // length Ncol
 "	__local "+ typeString+ " *Bcache,\n"
 "const int Nrow, const int Ncol, const int NpadC, const int NpadA,\n" 
 "const int NpadB, const int NlocalCache) {\n" +
 typeString + " Dout;\n"

"int Dglobal0 = get_global_id(0), Nglobal0 = get_global_size(0);\n"
"int Nglobal1 = get_global_size(1);\n"
"int Drow, Dcol, Dinner, DinnerStop, DrowNpadC, DrowNpadA;\n"
"int Dlocal0 = get_local_id(0);\n"
"int Nlocal0 = get_local_size(0);\n"
"int Dlocal1 = get_local_id(1);\n"
"int Nlocal1 = get_local_size(1);\n"
"int doCacheA = (get_local_id(1) == 0);\n"

// work items shoudl be Nglobal1 by Ncol, local 1 by Ncol

// cache first rows of B 
"for(Drow = Dlocal0; Drow < NlocalCache; Drow += Nlocal0){\n"
  "for(Dcol = Dlocal1; Dcol < Ncol; Dcol += Nlocal1){\n"
    "Bcache[Dcol + Drow*Ncol] = B[Dcol + Drow * NpadB];\n"
  "}\n"
"}\n"


// looped through rows which are all cached
"DinnerStop = min(Nrow, NlocalCache);\n"
"for(Drow = Dglobal0; Drow < DinnerStop; Drow+=Nglobal0){\n"
  "DrowNpadA= Drow * NpadA;\n"
  "DrowNpadC= Drow * NpadC;\n"
  "for(Dcol = get_global_id(1); Dcol < Ncol; Dcol += Nglobal1){\n"
    "Dout = 0.0;\n"
    "for(Dinner = 0; Dinner <= Drow; Dinner++){\n"
	  "if(doCacheA) {\n"
	    "Acache[Dlocal0] = A[Dinner + DrowNpadA];\n"
	  "}\n"
	  "barrier(CLK_LOCAL_MEM_FENCE);\n"
	  "Dout += Acache[Dlocal0] * Bcache[Dcol + Dinner * Ncol];\n"
	"}\n" // Dinner
  	"C[Dcol + DrowNpadC] = Dout;\n"
  "}\n" // Dcol
"}\n" //Drow
// rows which are not all cached
"for( ; Drow < Nrow; Drow+=Nglobal0){\n"
  "DrowNpadA= Drow * NpadA;\n"
  "DrowNpadC= Drow * NpadC;\n"
  "for(Dcol = get_global_id(1); Dcol < Ncol; Dcol += Nglobal1){\n"
    "Dout = 0.0;\n"
    // cached rows
    "for(Dinner = 0; Dinner < NlocalCache; Dinner++){\n"
	  "if(doCacheA) {\n"
	    "Acache[Dlocal0] = A[Dinner + DrowNpadA];\n"
	  "}\n"
	  "barrier(CLK_LOCAL_MEM_FENCE);\n"
	  "Dout += Acache[Dlocal0] * Bcache[Dcol + Dinner * Ncol];\n"
	"}\n" // Dinner
	// un-cached rows
    "for( ; Dinner <= Drow; Dinner++){\n"
	  "if(doCacheA) {\n"
	    "Acache[Dlocal0] = A[Dinner + DrowNpadA];\n"
	  "}\n"
	  "barrier(CLK_LOCAL_MEM_FENCE);\n"
	  "Dout += Acache[Dlocal0] * B[Dcol + Dinner * NpadB];\n"
	"}\n" // Dinner
	"C[Dcol + DrowNpadC] = Dout;\n"
  "}\n" // Dcol
 "}\n" //Drow
"}";

  return(result);
}





template <typename T> 
void multiplyLower(
	viennacl::matrix<T> &C,
	viennacl::matrix<T> &A,
	viennacl::matrix<T> &B,
	const int Nglobal0, const int Nlocal0, const int NlocalCache, const int ctx_id) {

	// the context
	viennacl::ocl::context ctx(viennacl::ocl::get_context(ctx_id));

	cl_device_type type_check = ctx.current_device().type();

	viennacl::ocl::program & my_prog = ctx.add_program(
		multiplyLowerTypeString<T>(),
		"my_kernel");

	viennacl::ocl::kernel 
		&multiplyLowerKernel = my_prog.get_kernel("multiplyLower");

	int Ncol = C.size2();
	multiplyLowerKernel.global_work_size(0, Nglobal0);
	multiplyLowerKernel.global_work_size(1, Ncol);
	multiplyLowerKernel.local_work_size(0, Nlocal0);
	multiplyLowerKernel.local_work_size(1, Ncol);

	viennacl::ocl::local_mem Acache(Nlocal0*sizeof(cl_double)); // T?
	viennacl::ocl::local_mem Bcache(NlocalCache*Ncol*sizeof(cl_double));


		// diagonals and diagTimesRowOfA
		viennacl::ocl::enqueue(multiplyLowerKernel(
			C, A, B,
			Acache, Bcache,
			(int) A.size1(), 
			Ncol, 
			(int) C.internal_size2(), 
			(int) A.internal_size2(), 
			(int) B.internal_size2(),
			(int) NlocalCache
			));
	
}


//' Multiply lower triangular matrix
//' 
//' Multiplies a lower triangular matrix by a rectangular matrix
//'
//' @param C output matrix
//' @param A lower triangular matrix
//' @param B rectangular matrix
//' @param Nglobal0 number of global work items
//' @param Nlocal0 number of local work items
//' @param NlocalCache rows of B to cache in local memory
//' @export
// [[Rcpp::export]]
SEXP multiplyLowerDouble(
	Rcpp::S4 C,
	Rcpp::S4 A,
	Rcpp::S4 B,
	int Nglobal0,
	int Nlocal0, 
	int NlocalCache) {

	const int ctx_id = INTEGER(C.slot(".context_index"))[0]-1;
	const bool BisVCL=1;

	std::shared_ptr<viennacl::matrix<double> > 
		AG = getVCLptr<double>(A.slot("address"), BisVCL, ctx_id);
	std::shared_ptr<viennacl::matrix<double> > 
		BG = getVCLptr<double>(B.slot("address"), BisVCL, ctx_id);
	std::shared_ptr<viennacl::matrix<double> > 
		CG = getVCLptr<double>(C.slot("address"), BisVCL, ctx_id);

	multiplyLower<double>(*CG, *AG, *BG, Nglobal0, Nlocal0, NlocalCache, ctx_id);	

	return Rcpp::wrap(0L);
}
