#include "geostatsgpu.hpp"
//#define DEBUG



template <typename T> 
std::string cholBatch(
  int N,
  int Npad,
  int NpadDiag,
  int Nmatrix,
  int NpadBetweenMatrices,
  int Ncache, // must exceed Nlocal0 * Nlocal1
  int Nlocal0, int Nlocal1, int Nlocal2, 
  int NlocalGroups0,int NlocalGroups1, int NlocalGroups2
) {

  std::string typeString = openclTypeString<T>();
  std::string result = "";

  if(typeString == "double") {
  	result += "\n#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n\n";
  }

  result +=
  "\n#define N " + std::to_string(N) + "//dimension of matrix\n"
  "#define Npad " + std::to_string(Npad) + "//internal number of columns\n"
  "#define NpadDiag " + std::to_string(NpadDiag) + "//internal columns for matrix holding diagonals\n"
  "#define Nmatrix " + std::to_string(Nmatrix) + "// number of matrices\n"
  "#define Ncache " + std::to_string(Ncache) + "//elements in internal cache\n"
  "#define NpadBetweenMatrices " + std::to_string(NpadBetweenMatrices) + "//extra rows between stacked matrices\n\n";


  result = result + 
"__kernel void cholGpu("
"	__global " + typeString + " *A," 
"	__global " + typeString + " *diag,"
"	const int Nstop){\n"

"// first dimension is rows"
"// second dimension is work items doing the same row"
"// third dimension matrix"


" const int doDiag = (get_global_id(0) == 0);// diagonals only done with first row\n" 
" const int doLocalDiagSum = (doDiag & get_local_id(1) == 0);\n"
" const int doGlobalDiagSum = (doDiag & get_global_id(1) == 0);\n"
"const int diagLocalAddIndex = get_local_id(1)*get_local_size(1)+get_local_id(2);\n"
"const int diagGlobalAddIndex = get_group_id(1)*get_num_groups(1)+get_group_id(2);\n"
"const int localAddIndex = get_local_id(0)*get_local_size(0)+get_local_id(2);\n"
"const int globalAddIndex = get_group_id(0)*get_group_size(0)+get_global_id(2);\n"

" __local " + typeString + " diagLocal[Ncache];//local cache of diagonals\n"
" __local " + typeString + " toAddLocal["+std::toString(
  max(Nlocal1*Nlocal2, // for diagonals
    Nlocal0*Nlocal2) // for off-diagonals
  )+"];//Nlocal0 * Nlcoal1\n"
" __global " + typeString + " toAddGlobal["+std::toString(
  max(
    NlocalGroups1,
    NlocalGroups0 * NlocalGroups2) + 
  "];//NlocalGroups1\n"


"	int Dcol, DcolNpad, Dcolm1;\n"
"	int Drow, Dk, DrowDk, Dmatrix, toAddHere;\n"
typeString + " DL, Ddiag, *AHere, *diagHere, *toAddGlobalHere;\n"

"for(Dcol =  0; Dcol < N; Dcol++) {"
"   DcolNpad = Dcol*Npad;\n"
"   Dcolm1 = Dcol - 1;"


// diagonal entry, use only first group
"if(get_global_id(0)==0){"

"for(Dmatrix = get_local_id(2);Dmatrix < Nmatrix;Dmatrix+= get_local_size(2)){"
"diagLocal[diagAddIndex]=0.0;"
"AHere = &A[Dmatrix*NpadBetweenMatrices];\n"
"diagHere = &diag[Dmatrix*NpadDiag];\n"
" for(Dk = get_global_id(1); Dk < Dcol; Dk += get_global_size(1)) {\n"
"   DL = AHere[DcolNpad + Dk];"
"   DL *= DL;"
"   diagLocal[diagLocalAddIndex] += diagHere[Dk] * DL;"
" }\n" // Dk

// local reduction
" barrier(CLK_LOCAL_MEM_FENCE);\n"
"if(get_local_id(1)==0){"
" for(Dk = 1; Dk < get_local_size(1); Dk++) {\n"
" diagLocal[diagLocalAddIndex] +=  diagLocal[diagLocalAddIndex + Dk];\n"
"}\n" //Dk
"}\n" //diagSum
" toAddGlobal[diagGlobalAddIndex] =  diagLocal[diagLocalAddIndex];\n"

// final reduction
" barrier(CLK_GLOBAL_MEM_FENCE);\n"
"if(get_global_id(1) == 0){"
"toAddGlobalHere = &toAddGlobal[Dmatrix * Nmatrix];"
" DL = toAddGlobalHere[0];"
" for(Dk = 1; Dk < get_group_size(1); Dk++) {\n"
"DL += toAddGlobalHere[Dk];"
"}\n" //Dk global diag sum
"   diagHere[Dcol] = AHere[Dcol + DcolNpad] - DL;\n"
"}\n" //diagGlobalSum

"#ifdef diagToOne\n"
"AHere[Dcol + DcolNpad] = 1.0;\n"
"#endif\n"

" barrier(CLK_GLOBAL_MEM_FENCE);\n"
"}\n" //Dmatrix

"	}\n" // doDiag

// off-diagonals
// TO DO: cache diag * A[Dcol, ]
"	for(Drow = Dcol+get_global_id(0); Drow < N; Drow += get_global_size(0)) {\n"

"for(Dmatrix = get_local_id(2);Dmatrix < Nmatrix;Dmatrix+= get_local_size(2)){"
"  AHere = &A[Dmatrix*NpadBetweenMatrices];\n"
"  diagHere = &diag[Dmatrix*NpadDiag];\n"
"	 DL = 0.0;\n"

"	 for(Dk = get_global_id(1); Dk < Dcolm1; Dk+=get_global_size(1)) {\n"
//			DL -= A[DrowDk] * diagLocal[Dk];"
"    DL += AHere[Dcol *  Npad + Dk] * AHere[Drow * Npad + Dk] * diagHere[Dk];"
			// DL -= A[Drow + Dk * Npad] * A[Dcol + DkNpad] * diag[Dk];"
"  }\n" // Dk
"  toAddLocal[addLocalIndex] = DL;\n"

// local reduction
"  barrier(CLK_LOCAL_MEM_FENCE);\n"
"  if(get_local_id(1) == 0){"
"    for(Dk = get_local_id(1)+1; Dk < get_local_size(1); Dk++) {\n"
"      toAddLocal[localAddIndex] +=  toAddLocal[localAddIndex + Dk];\n"
"    }\n" //Dk
"    toAddGlobal[diagGlobalAddIndex] =  toAddLocal[localAddIndex];\n"
"  }\n" //diagSum

// global reduction
"  barrier(CLK_GLOBAL_MEM_FENCE);\n"
"  if(get_global_id(1) == 0){\n"
"    toAddGlobalHere = &toAddGlobal[diagGlobalAddIndex];\n"
"    DL = diagLocal[addLocalIndex];\n"
"    for(Dk = 1; Dk < get_group_size(1); Dk++) {\n"
"       DL += toAddGlobalHere[Dk];\n"
"    }\n" //Dk global diag sum

"    AHere[Drow + DcolNpad] -= DL\n"
"    AHere[Drow + DcolNpad] /= diagHere[Dcol];\n"

"#ifdef upperToZero\n"
"    AHere[Drow*Npad + Dol] =0.0;\n"
"#endif\n"

"  }\n" //global Sum
"  barrier(CLK_GLOBAL_MEM_FENCE);\n"

"}\n" // Dmatrix loop
"}\n" // Drow
"}\n" // Dcol loop

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
