#include "geostatsgpu.hpp"
#define DEBUG

/* TO DO
- v1: one matrix per work group

    - for Dcol small
    - d0 matrix rows d1 multiple items on same row. local memory for reductions
    - cache A[Dcol, 1:Dcol] D[1:Dcol] in local memory
    - ... with a global cache for any overflow

- v2: small number of work items per matrix, final reductions on gpu

    - for N-Dcol large
    - work items rows by cols by matrix
    - loop from Dcol to Dcol + K
    - need K by (Dcol+K)/Ngroups[1] by 2 local cache for matrix multiplication
    - or K by Nlocal1 by (2 + c) local cache for matrix multiplication, cache c row groups locally, remainder global
    - need Nmatrix by K by (N-Dcol) by Ngroups global memory for final reduction
*/


template <typename T> 
std::string cholBatchKernelString( // V1
  int colStart,
  int colEnd,
  int N,
  int Npad,
  int NpadDiag,
  int Nmatrix,
  int NpadBetweenMatrices,
  int Ncache, 
  std::vector<int> Nlocal, // length 2
  bool allowOverflow
) {

  std::string typeString = openclTypeString<T>();
  std::string result = "";


  if(typeString == "double") {
  	result += "\n#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n\n";
  }

  result +=
  "\n//dimension of matrix\n#define N " + std::to_string(N) + "\n"
  "//column to start at\n#define colStart " + std::to_string(colStart) + "\n"
  "#define colEnd " + std::to_string(colEnd) + "\n"
  "//internal number of columns\n#define Npad " + std::to_string(Npad) + "\n"
  "//internal columns for matrix holding diagonals\n#define NpadDiag " + std::to_string(NpadDiag) + "\n"
  "#define Nmatrix " + std::to_string(Nmatrix) + "\n"
  "//elements in internal cache\n#define Ncache " + std::to_string(Ncache) + "\n"
  "//extra rows between stacked matrices\n#define NpadBetweenMatrices " + std::to_string(NpadBetweenMatrices) + "\n"
  "#define maxLocalItems " + std::to_string(Nlocal[0]*Nlocal[1]) + "\n\n";

if(allowOverflow) {
  result += "\n#define allowOverflow\n\n";
} else {
  result += "\n#define minDcolNcache Dcol\n" 
  "#define minDcolm1Ncache Dcolm1\n\n";
}
  
result += "__kernel void cholBatch(\n"
"	__global " + typeString + " *A,\n" 
"	__global " + typeString + " *diag\n"
"){\n"

//"const int colStart=0, colStop=;\n"
"const int localIndex = get_local_id(0)*get_local_size(1) + get_local_id(1);\n"
"const int NlocalTotal = get_local_size(0)*get_local_size(1);\n"
"const int localAddIndex = get_local_id(0)*get_local_size(1)+get_local_id(1);\n"

" __local " + typeString + " diagLocal[Ncache];//local cache of diagonals\n"
" __local " + typeString + " toAddLocal[maxLocalItems];\n"

"	int Dcol, DcolNpad, Dcolm1;\n"
"	int Drow, Dk, Dmatrix;\n";

if(allowOverflow) {
result += "	int minDcolm1Ncache, minDcolNcache;\n";
}

result +=  typeString + " DL;\n" 
  "__local " +  typeString + " diagDcol;\n" 
  "__global " + typeString + " *AHere, *AHereDcol, *AHereDrow, *diagHere;\n"


"for(Dmatrix = get_group_id(0); Dmatrix < Nmatrix; Dmatrix+= get_num_groups(0)){\n"

"diagHere = &diag[Dmatrix*NpadDiag];\n"
"AHere = &A[Dmatrix*NpadBetweenMatrices];\n"


"for(Dcol = colStart; Dcol < colEnd; Dcol++) {\n"

"  DcolNpad = Dcol*Npad;\n"
"  Dcolm1 = Dcol - 1;\n"
"  AHereDcol = &AHere[DcolNpad];\n"

"\n#ifdef allowOverflow\n"
"  minDcolNcache = min(Dcol, Ncache);\n"
"  minDcolm1Ncache = min(Dcolm1, Ncache);\n"
"#endif\n\n"

// diagonals
"  toAddLocal[localIndex]=0.0;\n"
"  for(Dk = localIndex; Dk < minDcolNcache; Dk += NlocalTotal) {\n"
"    DL = AHereDcol[Dk];\n"
"    diagLocal[Dk] = diagHere[Dk] * DL;\n"// cached A[Dcol, 1:Dcol] D[1:Dcol]
"    toAddLocal[localIndex] += diagLocal[Dk] * DL;\n"
"  }\n\n" // Dk
"\n#ifdef allowOverflow\n"
"  for(; Dk < Ncache; Dk += NlocalTotal) {\n"
"    DL = AHereDcol[Dk];\n"
"    toAddLocal[localIndex] += diagHere[Dk] * DL * DL;\n"
"  }\n\n" // Dk
"\n#endif\n"


// reduction on dimension 1
"barrier(CLK_LOCAL_MEM_FENCE);\n"
"if(get_local_id(1) == 0){"
" for(Dk = 1; Dk < get_local_size(1); Dk++) {\n"
" toAddLocal[localIndex] +=  toAddLocal[localIndex + Dk];\n"
"}\n" //Dk
"}\n" //get_local_id(1) == 0

// final reduction on dimension 0
"barrier(CLK_LOCAL_MEM_FENCE);\n"
"if(localIndex==0){\n"
"  for(Dk = localIndex + get_local_size(1); Dk < NlocalTotal; Dk+= get_local_size(1)) {\n"
"    toAddLocal[localIndex] +=  toAddLocal[Dk];\n"
"  }\n" //Dk

"  diagDcol = AHereDcol[Dcol] - toAddLocal[localIndex];\n"
"  diagHere[Dcol] = diagDcol;\n"
#ifdef DEBUG
"AHere[Dcol] = -toAddLocal[localIndex];\n"
#endif
"\n#ifdef diagToOne\n"
"AHere[Dcol] = 1.0;\n"
"#endif\n"
"}\n" //localIndex==0
"  barrier(CLK_LOCAL_MEM_FENCE);\n"
"  diagDcol = diagHere[Dcol];\n"

// off diagonals

"	for(Drow = Dcol+get_local_id(0)+1; Drow < N; Drow += get_local_size(0)) {\n"

"  AHereDrow = &AHere[Drow*Npad];\n"
"	 DL = 0.0;\n"

"	 for(Dk = get_local_id(1); Dk < minDcolm1Ncache; Dk+=get_local_size(1)) {\n"
"    DL += AHereDrow[Dk] * diagLocal[Dk];\n"
			// DL -= A[Drow + Dk * Npad] * A[Dcol + DkNpad] * diag[Dk];"
"  }\n" // Dk
"\n#ifdef allowOverflow\n"
"	 for(; Dk < Dcolm1; Dk+=get_local_size(1)) {\n"
"    DL += AHereDrow[Dk] * diagHere[Dk] * AHereDcol[Dk];\n"
"  }\n" // Dk
"\n#endif\n"
"  toAddLocal[localIndex] = DL;\n"

// local reduction
"barrier(CLK_LOCAL_MEM_FENCE);\n"
"if(get_local_id(1) == 0){"
"  DL = toAddLocal[localIndex];\n"
"  for(Dk = 1; Dk < get_local_size(1); Dk++) {\n"
"    DL +=  toAddLocal[localIndex + Dk];\n"
"  }\n" //Dk
"  AHereDrow[Dcol] = (AHereDrow[Dcol] - DL)/diagDcol;\n"
#ifdef DEBUG
"  AHereDcol[Drow] = - DL;\n" 
#endif
"}//get_local_id(1) == 0\n" 
"}//Drow\n" 

"  barrier(CLK_GLOBAL_MEM_FENCE);\n"
"} // Dcol loop\n"
"} // Dmatrix loop\n"
"}\n";
return(result);
}

template <typename T> 
int cholBatchVcl(
  viennacl::matrix<T> &A,
  viennacl::matrix<T> &D,
  const std::vector<int> &Nglobal,
  const std::vector<int> &Nlocal, 
  const int NlocalCache,
  const int ctx_id) {

  std::string cholClString = cholBatchKernelString<T>(
  0L, // start
//  A.size2(), // end
  2,
  A.size2(), // N
  A.internal_size2(), // Npad
  D.internal_size2(),
//  D.size1(), // Nmatrix
  2, 
  A.size2() * A.internal_size2(),// NpadBetweenMatrices,
  NlocalCache, 
  Nlocal,
  A.size2() > NlocalCache); // allow overflow

  viennacl::ocl::context ctx(viennacl::ocl::get_context(ctx_id));
  viennacl::ocl::program & my_prog = ctx.add_program(cholClString, "my_kernel");
  
#ifdef DEBUG
  
  Rcpp::Rcout << cholClString << "\n\n";
  
#endif  
  
  viennacl::ocl::kernel & cholKernel = my_prog.get_kernel("cholBatch");
  
  if(Nlocal[1] != Nglobal[1]) {
    Rf_warning("local and global work sizes should be identical for dimension 2");
  }

// dimension 0 is cell, dimension 1 is matrix
  cholKernel.global_work_size(0, (cl_uint) (Nglobal[0] ) );
  cholKernel.global_work_size(1, (cl_uint) (Nglobal[1] ) );

  cholKernel.local_work_size(0, (cl_uint) (Nlocal[0]));
  cholKernel.local_work_size(1, (cl_uint) (Nlocal[1]));

  // with opencl 2.0 can put this as a program-level variable
//  int Ngroups1 = ceil(Nglobal[1]/Nlocal[1]);
//  viennacl::matrix<T> finalReduction(Nglobal[0]*Nglobal[2]*Ngroups1);//size Nglobal[0]*Ngroups1*Nglobal[2]

  viennacl::ocl::enqueue(cholKernel(A, D));//,  2L));//A.size1() ));

  return 0L;
}


template<typename T> void cholBatchTemplated(
  Rcpp::S4 A,
  Rcpp::S4 D,
  const std::vector<int> &Nglobal,
  const std::vector<int> &Nlocal, 
  const int NlocalCache
) {

 
  // data
  const bool BisVCL=1;
  const int ctx_id = INTEGER(A.slot(".context_index"))[0]-1;
  std::shared_ptr<viennacl::matrix<T> > vclA = getVCLptr<T>(A.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > vclD = getVCLptr<T>(D.slot("address"), BisVCL, ctx_id);

  cholBatchVcl<T>(
    *vclA, *vclD, 
    Nglobal, 
    Nlocal,
    NlocalCache,
    ctx_id);

}


//[[Rcpp::export]]
void cholBatchBackend(
  Rcpp::S4 A,
  Rcpp::S4 D,
  std::vector<int> Nglobal,
  std::vector<int> Nlocal,
  int NlocalCache=500L) {


    Rcpp::traits::input_parameter< std::string >::type classVarR(RCPP_GET_CLASS(A));
    std::string precision_type = (std::string) classVarR;

  if(precision_type == "fvclMatrix") {
  cholBatchTemplated<float>(
    A, D, Nglobal, Nlocal,
    NlocalCache);
  } else if (precision_type == "dvclMatrix") {
  cholBatchTemplated<double>(
    A, D, Nglobal, Nlocal,
    NlocalCache);
  } else {
    Rcpp::warning("class of A must be fvclMatrix or dvclMatrix");
 }
}
