#include "geostatsgpu.hpp"
//#define DEBUG



template <typename T> 
std::string cholBatchKernelString( // V1
  int colStart,
  int colEnd,
  int N,
  int Npad,
  int NpadDiag,
  int Nmatrix,
  int NpadBetweenMatrices,
  std::vector<int> Ncache, 
  std::vector<int> Nlocal, // length 2
  bool allowOverflow
) {

  std::string typeString = openclTypeString<T>();
  std::string result = "";


  if(typeString == "double") {
  	result += "\n#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n\n";
  }

  result +=
  "\n#define N " + std::to_string(N) + "  //dimension of matrix\n"
  "#define colStart " + std::to_string(colStart) + "  //column to start at\n"
  "#define colEnd " + std::to_string(colEnd) + "\n"
  "//internal number of columns\n#define Npad " + std::to_string(Npad) + "\n"
  "//internal columns for matrix holding diagonals\n#define NpadDiag " + std::to_string(NpadDiag) + "\n"
  "#define Nmatrix " + std::to_string(Nmatrix) + "\n"
  "//elements in internal cache\n#define Ncache " + std::to_string(Ncache[0]) + "\n"
  "//extra rows between stacked matrices\n#define NpadBetweenMatrices " + std::to_string(NpadBetweenMatrices) + "\n"
  "#define maxLocalItems " + std::to_string(Nlocal[0]*Nlocal[1]) + "\n";


result += "\n__kernel void cholBatch(\n"
"	__global " + typeString + " *A,\n" 
"	__global " + typeString + " *diag\n"
"){\n"

" const int localIndex = get_local_id(0)*get_local_size(1) + get_local_id(1);\n"
" const int NlocalTotal = get_local_size(0)*get_local_size(1);\n"

" local " + typeString + " diagLocal[Ncache];//local cache of diagonals\n"

" local " + typeString + " toAddLocal[maxLocalItems];\n"

"	int Dcol, DcolNpad;\n"
"	int Drow, Dk, Dmatrix;\n";

if(allowOverflow) {
result += "	int minDcolNcache;\n";
}

result +=  typeString + " DL;\n" 
"  local " +  typeString + " diagDcol;\n" 
"  int AHere, AHereDcol, AHereDrow, diagHere;\n"
  

"for(Dmatrix = get_group_id(0); Dmatrix < Nmatrix; Dmatrix+= get_num_groups(0)){\n"

"diagHere = Dmatrix*NpadDiag;\n"
"AHere = Dmatrix*NpadBetweenMatrices;\n"

"for(Dcol = colStart; Dcol < colEnd; Dcol++) {\n"
"  DcolNpad = Dcol*Npad;\n"
"  AHereDcol = AHere+DcolNpad;\n"
//"  Dcolm1 = Dcol - 1;\n"
"  toAddLocal[localIndex]=0.0;\n";


if(allowOverflow) {
  result +="  minDcolNcache = min(Dcol, Ncache);\n"
  "  for(Dk = localIndex; Dk < minDcolNcache; Dk += NlocalTotal) {\n";
} else {
  result +="  for(Dk = localIndex; Dk < Dcol; Dk += NlocalTotal) {\n";
}


// diagonals
result += 
"    DL = A[AHereDcol+Dk];\n"
"    diagLocal[Dk] = diag[diagHere+Dk] * DL;\n"// cached A[Dcol, 1:Dcol] D[1:Dcol]
"    toAddLocal[localIndex] += diagLocal[Dk] * DL;\n"
"  }// Dk\n"; 


if(allowOverflow) {
  result +=
    "  for(Dk=minDcolNcache+localIndex; Dk < N; Dk += NlocalTotal) {\n"
"    DL = A[AHereDcol+Dk];\n"
"    toAddLocal[localIndex] += diag[diagHere+Dk] * DL * DL;\n"
"  }// Dk\n";
}


// reduction on dimension 1
result += "barrier(CLK_LOCAL_MEM_FENCE);\n"
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

"  diagDcol = A[AHereDcol+Dcol] - toAddLocal[localIndex];\n"
"  diag[diagHere+Dcol] = diagDcol;\n"
#ifdef DEBUG
"A[AHere+Dcol] = -toAddLocal[localIndex];\n"
//"AHereDcol[Dcol] = 10*Dcol;\n"
#endif
"\n#ifdef diagToOne\n"
"A[AHereDcol+Dcol] = 1.0;\n"
"#endif\n"
"}\n" //localIndex==0
"  barrier(CLK_LOCAL_MEM_FENCE);\n";
//"  diagDcol = diagHere[Dcol];\n"

// off diagonals
result +=
"	for(Drow = Dcol+get_local_id(0)+1; Drow < N; Drow += get_local_size(0)) {\n"

"  AHereDrow = AHere+Drow*Npad;\n"
"	 DL = 0.0;\n";


if(allowOverflow) {
  result +=
    "	 for(Dk = get_local_id(1); Dk < minDcolNcache; Dk+=get_local_size(1)) {\n";
} else {
  result +=
    "	 for(Dk = get_local_id(1); Dk < Dcol; Dk+=get_local_size(1)) {\n";
}

result +=
  "    DL += A[AHereDrow+Dk] * diagLocal[Dk];\n"
			// DL -= A[Drow + Dk * Npad] * A[Dcol + DkNpad] * diag[Dk];"
"  } // Dk\n";

if(allowOverflow) {
  result +=
    "	 for(minDcolNcache + get_local_id(1); Dk < Dcol; Dk+=get_local_size(1)) {\n"
"    DL += A[AHereDrow+Dk] * diag[diagHere+Dk] * A[AHereDcol+Dk];\n"
"  }// Dk\n"; 
}
result +=
  "  toAddLocal[localIndex] = DL;\n";

// local reduction
result +=
  "barrier(CLK_LOCAL_MEM_FENCE);\n"
"if(get_local_id(1) == 0){"
"  DL = toAddLocal[localIndex];\n"
"  for(Dk = 1; Dk < get_local_size(1); Dk++) {\n"
"    DL +=  toAddLocal[localIndex + Dk];\n"
"  }\n" //Dk
"  A[AHereDrow+Dcol] = (A[AHereDrow+Dcol] - DL)/diagDcol;\n"
#ifdef DEBUG
"  A[AHereDcol+Drow] = DL;\n" 
#endif
"}//get_local_id(1) == 0\n\n"; 

result +=
  "}//Drow\n"
  "  barrier(CLK_GLOBAL_MEM_FENCE);\n"
"} // Dcol loop\n"
"} // Dmatrix loop\n\n"

"}\n";
return(result);
}

template <typename T> 
int cholBatchVcl(
  viennacl::matrix<T> &A,
  viennacl::matrix<T> &D,
  const std::vector<int> &Nglobal,
  const std::vector<int> &Nlocal, 
  const std::vector<int> &NlocalCache,
  const int ctx_id) {

  std::string cholClString = cholBatchKernelString<T>(
  0L, // start
  A.size2(), // end
  A.size2(), // N
  A.internal_size2(), // Npad
  D.internal_size2(),
  D.size1(), // Nmatrix
  A.size2() * A.internal_size2(),// NpadBetweenMatrices,
  NlocalCache, 
  Nlocal,
  A.size2() > NlocalCache[0]); // allow overflow

  viennacl::ocl::context ctx(viennacl::ocl::get_context(ctx_id));
  viennacl::ocl::program & my_prog = ctx.add_program(cholClString, "my_kernel");
  
#ifdef DEBUG
  
  Rcpp::Rcout << cholClString << "\n\n";
  
#endif  
  
  viennacl::ocl::kernel & cholKernel = my_prog.get_kernel("cholBatch");
  
  if(Nlocal[1] != Nglobal[1]) {
    Rf_warning("local and global work sizes should be identical for dimension 2, ignoring global");
  }

// dimension 0 is cell, dimension 1 is matrix
  cholKernel.global_work_size(0, (cl_uint) (Nglobal[0] ) );
  cholKernel.global_work_size(1, (cl_uint) (Nlocal[1] ) );

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
  const std::vector<int> &NlocalCache
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
  std::vector<int> NlocalCache) {


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
