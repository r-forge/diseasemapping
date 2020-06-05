#include "gpuRandom.hpp"
//#define DEBUG

// C = A^T A or A^T D A or A^T D^(-1) A 
// TO DO: iterate cache


template <typename T> 
std::string crossprodBatchString(
    const int Nrow, 
    const int Ncol,
    const int Nmatrix, 
    const int NpadC, 
    const int NpadA,
    const int NpadD, // set to zero to omit D
    const int invertD, // set to 1 for A^T D^(-1) A
    const int NpadBetweenMatricesC,
    const int NpadBetweenMatricesA,
    const int NlocalCacheA, // numbers of rows to cache of A
    const std::vector<int> Nlocal// cache a Nlocal[0] by Nlocal[1] submatrix of C
  ) { 
    
 /*
  * global groups col by matrix
  * local items inner by row
  * 
  */
 
  std::string typeString = openclTypeString<T>();
  std::string result = "";
  
  if(typeString == "double") {
    result += "\n#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n\n";
  }
  result += 
    "\n#define Nrow " + std::to_string(Nrow) + "\n"    
    "#define Ncol " + std::to_string(Ncol) + "\n"    
    "#define Nmatrix " + std::to_string(Nmatrix) + "\n"    
    "#define NpadC " + std::to_string(NpadC) + "\n"    
    "#define NpadA " + std::to_string(NpadA) + "\n"    
    "#define NpadD " + std::to_string(NpadD) + "\n"    
    "#define NpadLocal " + std::to_string(Nlocal[1]) + "\n"    
    "#define NpadBetweenMatricesC " + std::to_string(NpadBetweenMatricesC) + "\n"    
    "#define NpadBetweenMatricesA " + std::to_string(NpadBetweenMatricesA) + "\n"    
    "#define NrowStop " + std::to_string(std::min(NlocalCacheA, Nrow)) + "\n"    
    "#define NlocalCacheA "  + std::to_string(NlocalCacheA) + "\n\n";   
    
    result +=
    "__kernel void crossprodBatch(\n"
    "	__global " + typeString+ " *C,\n"
    "	__global "+ typeString+ " *B) {\n\n"
      
    "local " + typeString + " Acache[" + 
      std::to_string(NlocalCacheA) + "];\n" 
    "local " + typeString + " Ccache[" + 
      std::to_string(Nlocal[0] * Nlocal[1]) + "];\n" +

  typeString + "Cout, Ctemp;\n"
  "event_t ev\n;"
  "int AHere, CHere;\n"
  "int Dmatrix, Drow, Dcol, DrowNpadC, DcolNpadA, DrowNpadA, Dinner;\n"
  "const int AHereInc = get_num_groups(1)*NpadBetweenMatricesA;\n"
  "const int CHereInc = get_num_groups(1)*NpadBetweenMatricesC;\n"
  "const int DlocalInc = get_local_size(0)*NpadLocal;\n"
  "const int DlocalIncDiag = DlocalInc * DrowNpadA;\n"
  "const int DcolNpadAInc = get_num_groups(0)*NpadA;\n"
  "const int DrowNpadCInc = get_local_size(1)*NpadC;\n"
  "const int DrowNpadAInc = get_local_size(1)*NpadA;\n"
  "const int doLocalSum = (get_local_id(0)==0);\n"
  "const int doFinalSum = (get_local_id(0)==0 & get_local_id(1)==0);\n"
  "const int cacheIndex = get_local_id(1)+NpadLocal*get_local_id(0);\n";
  
  if(NpadD) {
    result += "int DHere;\n"
    "int DHereInc = get_num_groups(1)*NpadD;\n";
  }
  
      result +=  "\n\n"
  "for(Dmatrix = get_group_id(1),\n"
  "    AHere = Dmatrix * NpadBetweenMatricesA,\n";
  if(NpadD) {
    "    DHere = Dmatrix * NpadD,\n";
  }
  result +=  
    "    CHere = Dmatrix * NpadBetweenMatricesC;\n"
  "  Dmatrix < Nmatrix;\n"
  "  Dmatrix += get_num_groups(1),\n"
  "    AHere += AHereInc,\n";
  if(NpadD) {
    "    DHere += DHereInc,\n";
  }
  result +=  
  "    CHere += CHereInc){\n";

  result +=  "\n"
  "  for(Dcol = get_group_id(0),\n"
  "      DcolNpadA = AHere + Dcol * NpadA;\n"
  "    Dcol < Ncol;\n"
  "    Dcol += get_num_groups(0),\n"
  "      DcolNpadA += DcolNpadAInc){\n";
  
  
  /*\
   * TO DO: 
   * switch Drow, Dinner loops, Nrounds cache
   */
  
  result +=  "\n"
  "    ev=async_work_group_strided_copy (Acache, &A[AHere],\n"
  "      NlocalCacheA, NpadA, 0);\n"
  "    wait_group_events (1, &ev);\n";

  // diagonal
  result += 
    "    Cout=0.0;\n"
    "    for(Dinner=cacheIndex;Dinner < NrowStop; Dinner += DlocalInc){\n";
  if(NpadD) {
    if(invertD) {
      result += 
        "      Ctemp = Acache[Dinner] / D[DHere+Dinner];\n";
    } else {
      result += 
        "      Ctemp = Acache[Dinner] * D[DHere+Dinner];\n";
    }
    result += 
      "      Cout += Acache[Dinner] * Ctemp;\n"
      "      Acache[Dinner] = Ctemp;\n";
  } else {
    result += 
      "      Cout += Acache[Dinner]*Acache[Dinner];\n";
  }
  result += 
    "    }\n"
    "    for(Dinner = NrowStop + cacheIndex,\n"
    "          DrowNpadA = DcolNpadA + Dinner * NpadA;\n"
    "        Dinner < Nrow;\n"
    "        Dinner += DlocalInc, DrowNpadA += DlocalIncDiag){\n"
    "      Ctemp = A[DrowNpadA];\n"
    "      Cout += Ctemp * Ctemp;\n"
    "    }\n"
    "    Ccache[cacheIndex] = Cout;\n"
    "    barrier(CLK_LOCAL_MEM_FENCE);\n"
    "    if(doFinalSum) {\n"
    "     for(Drow = 0,DrowNpadA=0; Drow < get_local_size(1);\n"
    "         Drow++,DrowNpadA += NpadLocal){\n"
    "       for(Dinner = 0;Dinner < get_local_size(0); Dinner++){\n"
    "         Cout += Ccache[DrowNpadA + Dinner];\n"
    "       }\n"
    "     }\n"
    "     C[DcolNpadC + Dcol] = Cout;\n"
    "    }\n"
    "    barrier(CLK_LOCAL_MEM_FENCE);\n";
    
    // off-diagonals
  result += 
  "    for(Drow = 1 + Dcol + get_local_id(1),\n"
  "        DrowNpadC = CHere + Drow * NpadC,\n"
  "        DrowNpadA = AHere + Drow * NpadA;\n"
  "      Drow < Ncolm1;\n"
  "      Drow += get_local_size(1),\n" 
  "        DcolNpadC += DrowNpadCInc,\n"
  "        DrowNpadA += DrowNpadAInc ){\n";
  
  result +=  "\n"
  "      Cout=0.0;\n"
  "      for(Dinner = get_local_id(0);\n"
  "        Dinner < Nstop;\n"
  "        Dinner += get_local_size(0)){\n";
  result += 
    "      Cout += A[DrowNpadA + Dinner] * Acache[Dinner]\n";
  
  "      for(Dinner < Nstop + get_local_id(0);\n"
  "        Dinner < Nrow;\n"
  "        Dinner += get_local_size(0)){\n";

  if(NpadD) {
    if(invertD) {
      result += 
        "      Cout += A[DcolNpadA + Dinner] * A[DrowNpadA + Dinner] / D[DHere+Dinner];\n";
    } else {
      result += 
        "      Cout += A[DcolNpadA + Dinner] * A[DrowNpadA + Dinner] * D[DHere+Dinner];\n";
    }
  } else {
    result += 
      "      Cout += A[DrowNpadA + Dinner] *\n"
      "        A[DcolNpadA + Dinner]\n";
  }  
    
    result += 
      "      }// Dinner\n";
    result +=       
      "      Ccache[cacheIndex] = Cout;\n"
      "      barrier(CLK_LOCAL_MEM_FENCE);\n"
      "      if(doLocalSum){"
      "for(Dinner = 1;Dinner < get_local_size(0);Dinner++){"
      "  Ccache[cacheIndex] += Ccache[cacheIndex + Dinner * NpadLocal;\n"
      "}\n"
      "C[DrowNpadC + Dcol] = Ccache[cacheIndex];\n"
      "      }\n"
      "      barrier(CLK_LOCAL_MEM_FENCE);\n";
      
    result += 
      "    }// Drow\n";
    result += 
      "  }// Dcol\n";
    result += 
      "}// Dmatrix\n";
    result += 
      "}// function";
  
  return(result);
}

template <typename T> 
void crossprodBatch(
    viennacl::matrix<T> &C,
    viennacl::matrix<T> &A,
    viennacl::matrix<T> &D,
    const int invertD,
    std::vector<int> Nglobal,
    std::vector<int> Nlocal,
    const int NlocalCache, 
    const int ctx_id) {
  
  
  const int Ncol = A.size2();
  const int Nmatrix = C.size1()/Ncol;
  const int Nrow = A.size1()/Nmatrix;

  
  
  // the context
  viennacl::ocl::context ctx(viennacl::ocl::get_context(ctx_id));
  
  cl_device_type type_check = ctx.current_device().type();
  
  const int NpadC, 
  const int NpadA,
  const int NpadD, // set to zero to omit D
  const int invertD, // set to 1 for A^T D^(-1) A
  const int NpadBetweenMatricesC,
  const int NpadBetweenMatricesA,
  
  std::string clString =  crossprodBatchString<T>(  
    Nrow, 
    Ncol, // ncol
    Nmatrix,
    C.internal_size2(), 
    A.internal_size2(), 
    D.internal_size2(),
    1, // A^T D^(-1) A
    C.internal_size2()*Nrow,//NpadBetweenMatricesC,
    A.internal_size2()*Nrow,//NpadBetweenMatricesA,
    NlocalCache,
    Nlocal);
  
#ifdef DEBUG
  
  Rcpp::Rcout << clString << "\n\n";
  
#endif  
  
  
  
  viennacl::ocl::program & my_prog = ctx.add_program(
    clString, "my_kernel");
  
  viennacl::ocl::kernel & multiplyKernel = my_prog.get_kernel("crossprodBatch");
  
  multiplyKernel.global_work_size(0, Nglobal[0]);
  multiplyKernel.global_work_size(1, Nglobal[1]);

  multiplyKernel.local_work_size(0, Nlocal[0]);
  multiplyKernel.local_work_size(1, Nlocal[1]);

  // diagonals and diagTimesRowOfA
  viennacl::ocl::enqueue(multiplyKernel(
      C, A, D, B));
  
}



template <typename T> 
SEXP crossprodBatchTyped(
    Rcpp::S4 CR,
    Rcpp::S4 AR,
    Rcpp::S4 DR,
    const int invertD,
    Rcpp::IntegerVector NglobalR,
    Rcpp::IntegerVector NlocalR, 
    const int NlocalCache) {
  
  std::vector<int> Nglobal = Rcpp::as<std::vector<int> >(NglobalR);
  std::vector<int> Nlocal = Rcpp::as<std::vector<int> >(NlocalR);
  
  const int ctx_id = INTEGER(CR.slot(".context_index"))[0]-1;
  const bool BisVCL=1;
  
  
  
  std::shared_ptr<viennacl::matrix<T> > 
    AG = getVCLptr<T>(AR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > 
    CG = getVCLptr<T>(CR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > 
    DG = getVCLptr<T>(DR.slot("address"), BisVCL, ctx_id);
  
  crossprodBatch<T>(*CG, *AG, *DG, invertD,
                                Nglobal, Nlocal, NlocalCache, ctx_id);
  
  return Rcpp::wrap(0L);
  
}


//' Multiply crossproduct matrices
//' 
//' Computes C = t(A) D A
//'
//' @param C output matrices, stacked row-wise
//' @param A rectangular matrices
//' @param D rectangular matrix, columns are diagonals
//' @param invertD set to 1 for C = t(A) D^(-1) A
//' @param Nglobal vector of number of global work items
//' @param Nlocal vector of number of local work items
//' @param NlocalCache elements in local cache
//' @export
// [[Rcpp::export]]
SEXP crossprodBatchBackend(
    Rcpp::S4 C,
    Rcpp::S4 A,
    Rcpp::S4 D,
    const int invertD,
    Rcpp::IntegerVector Nglobal,
    Rcpp::IntegerVector Nlocal, 
    const int NlocalCache) {
  
  SEXP result;
  
  Rcpp::traits::input_parameter< std::string >::type classVarR(RCPP_GET_CLASS(C));
  std::string precision_type = (std::string) classVarR;
  
  
  if(precision_type == "fvclMatrix") {
    result = crossprodBatchTyped<float>(C, A, D, invertD, Nglobal, Nlocal, NlocalCache);
  } else if (precision_type == "dvclMatrix") {
    result = crossprodBatchTyped<double>(C, A, D, invertD, Nglobal, Nlocal, NlocalCache);
  } else {
    result = Rcpp::wrap(1L);
  }
  return(result);
  
}

