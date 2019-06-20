#include "geostatsgpu.hpp"
#define DEBUG



template <typename T> 
std::string cholBatchKernelString(
  int N,
  int Npad,
  int NpadDiag,
  int Nmatrix,
  int NpadBetweenMatrices,
  int Ncache, // must exceed Nlocal0 * Nlocal1
  std::vector<int> Nlocal,
  std::vector<int> Nglobal
) {

  std::string typeString = openclTypeString<T>();
  std::string result = "";

//  int Ngroups0 = ceil(Nglobal[0]/Nlocal[0]);
  int Ngroups1 = ceil(Nglobal[1]/Nlocal[1]);


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
"\n __global " + typeString + " toAddGlobal["+std::to_string(Nglobal[0]*Ngroups1*Nglobal[2]) + "];//Nglobal0*Ngroups1*Nglobal2\n\n"
"__kernel void cholBatch(\n"
"	__global " + typeString + " *A,\n" 
"	__global " + typeString + " *diag,\n"
" const int colStart,\n"
"	const int colStop){\n"

"// first dimension is rows\n"
"// second dimension is work items doing the same row\n"
"// third dimension matrix\n"


"const int doDiag = get_global_id(2)==0 & get_group_id(1)==0;\n"
"const int localIndex = get_local_id(0)*get_local_size(1) + get_local_id(1);\n"
"const int NlocalTotal = get_local_size(0)*get_local_size(1);\n"
"const int localAddIndex = get_local_id(0)*get_local_size(1)+get_local_id(1);\n"
"const int globalAddIndex = get_global_size(0)*get_num_groups(1)*get_global_id(2) \n"
"    + get_global_id(0)*get_num_groups(1);\n"

" __local " + typeString + " diagLocal[Ncache];//local cache of diagonals\n"
" __local " + typeString + " toAddLocal["+std::to_string(Nlocal[0]*Nlocal[1])+"]; \n"


"	int Dcol, DcolNpad, Dcolm1;\n"
"	int Drow, Dk, Dmatrix, DmatrixWithPad;\n" +
typeString + " DL, diagDcol;\n" 
"__global " + typeString + " *AHere, *AHereDrow, *diagHere;\n"


"for(Dcol =  colStart; Dcol < colStop; Dcol++) {\n"
"   DcolNpad = Dcol*Npad;\n"
"   Dcolm1 = Dcol - 1;\n"


// diagonal entry, use only first group
"if(doDiag){\n"

"for(Dmatrix = get_group_id(0); Dmatrix < Nmatrix; Dmatrix+= get_num_groups(0)){\n"
"toAddLocal[localIndex]=0.0;\n"
"AHere = &A[Dmatrix*NpadBetweenMatrices + DcolNpad];\n"
"diagHere = &diag[Dmatrix*NpadDiag];\n"
"  for(Dk = localIndex; Dk < Dcol; Dk += NlocalTotal) {\n"
"    DL = AHere[Dk];\n"
"    DL *= DL;\n"
"    toAddLocal[localIndex] += diagHere[Dk] * DL;\n"
"  }\n\n" // Dk

// reduction on dimension 2
"barrier(CLK_LOCAL_MEM_FENCE);\n"
"if(get_local_id(1) == 0){"
" for(Dk = 1; Dk < get_local_size(1); Dk++) {\n"
" toAddLocal[localIndex] +=  toAddLocal[localIndex + Dk];\n"
"}\n" //Dk
"}\n" //diagSum

// final reduction on dimension 1
"barrier(CLK_LOCAL_MEM_FENCE);\n"
"if(localIndex==0){\n"
"  for(Dk = localIndex + get_local_size(1); Dk < NlocalTotal; Dk+= get_local_size(1)) {\n"
"    toAddLocal[localIndex] +=  toAddLocal[Dk];\n"
"  }\n" //Dk

"diagHere[Dcol] = AHere[Dcol] - toAddLocal[localIndex];\n"

"\n#ifdef diagToOne\n"
"AHere[Dcol] = 1.0;\n"
"#endif\n"


"}\n" //doFinalDiagSum

"barrier(CLK_LOCAL_MEM_FENCE);\n"

"}\n" //Dmatrix

"}\n" // doDiag
"barrier(CLK_GLOBAL_MEM_FENCE);\n"

// off-diagonals
"for(Dmatrix = get_global_id(2);Dmatrix < Nmatrix;Dmatrix+= get_global_size(2)){\n"
"  DmatrixWithPad = Dmatrix*NpadBetweenMatrices;\n"
"  AHere = &A[DmatrixWithPad + DcolNpad];\n"
"  diagHere = &diag[Dmatrix*NpadDiag];\n"
"  diagDcol = diagHere[Dcol];\n"

// TO DO: cache diag * A[Dcol, ]

"	for(Drow = Dcol+get_global_id(0); Drow < N; Drow += get_global_size(0)) {\n"

"  AHereDrow = &A[DmatrixWithPad + Drow*Npad];\n"
"	 DL = 0.0;\n"

"	 for(Dk = get_global_id(1); Dk < Dcolm1; Dk+=get_global_size(1)) {\n"
"    DL += AHereDrow[Dk] * AHere[Dk] * diagHere[Dk];\n"
			// DL -= A[Drow + Dk * Npad] * A[Dcol + DkNpad] * diag[Dk];"
"  }\n" // Dk
"  toAddLocal[localIndex] = DL;\n"

// local reduction
"  barrier(CLK_LOCAL_MEM_FENCE);\n"
"  if(get_local_id(1) == 0){\n"
"    for(Dk = 1; Dk < get_local_size(1); Dk++) {\n"
"      toAddLocal[localAddIndex] +=  toAddLocal[localAddIndex + Dk];\n"
"    }\n" //Dk
"  toAddGlobal[globalAddIndex] = toAddLocal[localIndex];\n"
"  }\n"  

// global reduction
"  barrier(CLK_GLOBAL_MEM_FENCE);\n"
"  if(get_global_id(1) == 0){\n"
"    toAddGlobalHere = &toAddGlobal[globalAddIndex];\n"
"    DL = toAddLocal[localIndex];\n"
"    for(Dk = 1; Dk < get_num_groups(1); Dk++) {\n"
"       DL += toAddGlobalHere[Dk];\n"
"    }\n" //Dk global diag sum

"    AHereDrow[Dcol] = (AHereDrow[Dcol] - DL)/diagDcol;\n"

"#ifdef upperToZero\n"
"    AHere[Drow] =0.0;\n"
"#endif\n"

"  }\n" //global reduction

"  barrier(CLK_GLOBAL_MEM_FENCE);\n"
"}\n" // Drow

"  barrier(CLK_GLOBAL_MEM_FENCE);\n"
"}\n" // Dmatrix loop

"  barrier(CLK_GLOBAL_MEM_FENCE);\n"
"}\n" // Dcol loop

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

  const int Nmatrix = D.size1();


  std::string cholClString = cholBatchKernelString<T>(
  A.size1(),
  A.internal_size2(),
  D.internal_size2(),
  Nmatrix,
  Nmatrix * A.internal_size2(),// NpadBetweenMatrices,
  NlocalCache, 
  Nlocal,
  Nglobal);

  viennacl::ocl::context ctx(viennacl::ocl::get_context(ctx_id));
  viennacl::ocl::program & my_prog = ctx.add_program(cholClString, "my_kernel");
  viennacl::ocl::kernel & cholKernel = my_prog.get_kernel("cholBatch");


// dimension 0 is cell, dimension 1 is matrix
  cholKernel.global_work_size(0, (cl_uint) (Nglobal[0] ) );
  cholKernel.global_work_size(1, (cl_uint) (Nglobal[1] ) );
  cholKernel.global_work_size(2, (cl_uint) (Nglobal[3] ) );

  cholKernel.local_work_size(0, (cl_uint) (Nlocal[0]));
  cholKernel.local_work_size(1, (cl_uint) (Nlocal[1]));
  cholKernel.local_work_size(2, (cl_uint) (Nlocal[2]));


#ifdef DEBUG

Rcpp::Rcout << cholClString << "\n\n";

#endif  

//  viennacl::ocl::enqueue(cholKernel(A, D, 0,  A.size1()));

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
  int NlocalCache=1000L) {


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
