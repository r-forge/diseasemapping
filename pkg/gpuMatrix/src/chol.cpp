
#include "gpuMatrix.hpp"



template <typename T> 
std::string cholMulti() {

  std::string typeString = openclTypeString<T>();
  std::string result = "";

  if(typeString == "double") {
  	result += "\n#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n\n";
  }
  result = result + 
"__kernel void cholGpu("
"	__global " + typeString + " *A," 
"	__global " + typeString + " *diag,"
"	__local " + typeString + " *diagLocal,"
"	const int N, "
"	const int Npad, const int NpadDiag,\n",
" const int NpadAA," //
"	const int Nstop, const int Nmatrix){\n"
// first dimension is rows, 

// second dimension is work items
// doing the same row

// third dimension matrix

"	const int Dglobal0 = get_global_id(0);" // row
" const int Nsize0 = get_global_size(0);"
"	const int Dglobal1 = get_global_id(1);" // group within row
" const int Nsize1 = get_global_size(1);"

" const int doDiag = (Dglobal0 == 0);" // diagonals only done with first row
" const int doDiagSum = (Dglobal0 == 0 & Dglobal1 == 0);\n"

// groups work on the same row of the same matrix
// only dimension 1 is different 
" const int Dlocal = get_local_id(1);\n"
" const int NlocalSize = get_local_size(1);\n"
" const int doDiagLocal = (get_local_id(0) == 0);\n"
// get_local_size(2) must be 1

"	int Dcol=0, Dcolm1;"
"	int DcolNpad = Dcol*Npad;"

"	__private int Drow, Dk, DrowDk, Dmatrix;"
"	__private " + typeString + " DL, Ddiag, *AHere, *diagHere;\n"

"for(Dmatrix =  get_global_id(2); Dmatrix < Nmatrix; Dmatrix += get_global_size(2)) {"

"AHere = &A[Dmatrix*NpadAA];\n"
"diagHere = &diag[Dmatrix*NpadDiag];\n"

"for(Dcol = 0; Dcol < Nstop; Dcol++){\n"
"   Dcolm1 = Dcol-1;\n"
"   DcolNpad = Dcol*Npad;\n"
// diagonal entry, use only first group
" if(doDiag){\n"
"	Ddiag = 0.0;\n"
"	for(Dk = Dglobal1; Dk < Dcol; Dk += Nsize1) {\n"
"		DL = AHere[Dcol + Dk * Npad];"
"		DL *= DL;"
"		diagLocal[Dlocal] += diagHere[Dk] * DL;"
"	}\n" // Dk
"	}\n" // doDiag
"	barrier(CLK_LOCAL_MEM_FENCE);\n"
// add 'em up to get diagonal
"	if(doDiagSum) {\n"
"		Ddiag = 0.0;\n"
"		for(Dk = 0; Dk < NlocalSize; Dk++) {\n"
"			Ddiag += diagLocal[Dk];\n"
"		}\n"
  // Ddiag is now diag[Dcol]
  // copy it to global memory
"		diagHere[Dcol] = AHere[Dcol + DcolNpad] - Ddiag;\n"
"   AHere[Dcol + DcolNpad] = 1.0;"
"	}\n" // doDiagSum
"	barrier(CLK_GLOBAL_MEM_FENCE);\n" // all items need the diagonal
" Ddiag = diagHere[Dcol];" 
 
	// diagLocal =  diag * A[Dcol, ]
" if(doDiagLocal){\n"
"	for(Dk = Dlocal; Dk < Dcol; Dk += NlocalSize) {\n"
"		DL = Ahere[Dcol + Dk * Npad];\n"
"		diagLocal[Dk] = diagHere[Dk] * DL;\n"
"	}\n"
" }\n" // doDiag
"	barrier(CLK_LOCAL_MEM_FENCE);\n" 


"	for(Drow = (N-Dglobal-1); Drow > Dcol; Drow -= Nsize0) {\n"

"		DL = AHere[Drow + DcolNpad];\n"

"		DrowDk = Dcolm1*Npad;\n"
"		for(Dk = Dcolm1; Dk >= 0; Dk--) {\n"
//			DL -= A[DrowDk] * diagLocal[Dk];"
"			DrowDk -= Npad;"
"			DL -= AHere[Drow + Dk * Npad] * diagLocal[Dk];"
			// DL -= A[Drow + Dk * Npad] * A[Dcol + DkNpad] * diag[Dk];"
"		}\n" // Dk
"		AHere[Drow + DcolNpad] = DL / Ddiag;"

"	}\n" // Drow

" barrier(CLK_GLOBAL_MEM_FENCE);\n" // all rows must be done before advancing to next column

"}\n" // Dcol loop
" barrier(CLK_GLOBAL_MEM_FENCE);\n" 
"}\n" // Dmatrix loop
"}\n";


template <typename T> 
int cholMulti(
  viennacl::matrix<T> *A,
  viennacl::vector<T> *D,
  const int Nmatrix,
  const IntegerVector Nglobal,
  const IntegerVector Nlocal, 
  const int NlocalCache,
  const int ctx_id) {

  const bool BisVCL=1;

  double **AA, **DD;
  

  return 0L;
}
