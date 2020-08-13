
#include "gpuRandom.hpp"
using namespace Rcpp;

/***
//#define DEBUG

// C = A B   C_i is M by N, A_i is M by K, B_i is K by N  (if no transposition)
// work items w1, w2, w3, i.e 32 by 32 by 4
// for now local items 1x1x1
// suppose matrices C are 120 by 140, 10 matrices
// item w1, w2, 0  will do matrices 0,4,8.  
// item w1, w2, 3 will do matrices 3, 7

// or C = A E B, if NpadE > 0

template <typename T> 
std::string gemmBatchString(
    const int M, // rows of Ai
    const int N, // columns of Bi
    const int K, // common size of Ai Bi
    const int NpadA, 
    const int NpadB, 
    const int NpadC,
    const int z, //number of batches//    const int Ntotal, //total columns of B
    const int NpadMatrixA, 
    const int NpadMatrixB, 
    const int NpadMatrixC,
    const int need_transpose) { 
  

  // Ncache must be bigger than local_size(1) and local_size(2)
  const int Ncache = 256;

  std::string typeString = openclTypeString<T>();
  std::string result = "";
  
  if(typeString == "double") {
    result += "\n#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n\n";
  }
  
  result += 
    "#define M " + std::to_string(M) + "\n" // rows 
    "#define N " + std::to_string(N) + "\n" // columns
    "#define K " + std::to_string(K) + "\n" // inter
    "#define z " + std::to_string(z) + "\n" // matrices
//    "#define Ntotal " + std::to_string(Ntotal) + "\n" // 
    "#define NpadA " + std::to_string(NpadA) + "\n"    
    "#define NpadB " + std::to_string(NpadB) + "\n"
    "#define NpadC " + std::to_string(NpadC) + "\n"
    "#define NpadMatrixA " + std::to_string(NpadMatrixA) + "\n"    
    "#define NpadMatrixB " + std::to_string(NpadMatrixB) + "\n"  
    "#define NpadMatrixC " + std::to_string(NpadMatrixC) + "\n"
    "#define Ncache " + std::to_string(Ncache) + "\n"; 
  
  result += "__kernel void gemm( __global "  + typeString+ "* A,\n"
                               " __global "  + typeString+ "* B,\n"
                               " __global " + typeString+ "* C){\n";
                              
    
  result += typeString + " acc;\n"
           "int Dmatrix, DmatrixA, DmatrixB;\n"
           "int Drow, Dcol, DrowC0, Dinner;\n";
  result += "int DrowLocal, DcolLocal, DinnerA, DinnerB;\n";
  result += "event_t wait;\n";

  result += "local " + typeString + " localCacheA[Ncache], localCacheB[Ncache];\n";



  result += 
  "for (Dmatrix=get_global_id(2);\n"
  "     Dmatrix < z;"
  "     Dmatrix += get_global_size(2)) {\n"
// for(DmatrixCol=0;DmatrixCol < NcolBatchA;DmatrixCol++){
  "  DmatrixA = Dmatrix * NpadMatrixA;\n"
  "  DmatrixB = Dmatrix * NpadMatrixB;\n";

 result += 
// DrowLocal is row for work item (0,0) 
  "  for (Drow = get_global_id(0), DrowLocal = get_group_id(0) * get_local_size(0);"
// "    Drow < M;"
  "    DrowLocal < M;\n" // to synchronise local memory, all work items must got through each cycle
  "    Drow += get_global_size(0),DrowLocal += get_global_size(0)) {\n"
  "    DrowC0 = Dmatrix * NpadMatrixC + Drow*NpadC;\n";
          
          
  result += 

// DcolLocal is column for work item (0,0) 
  "    for (Dcol = get_global_id(1), DcolLocal = get_group_id(1) * get_local_size(1);\n"
//  "      Dcol < Ntotal;"SHOULDNT THIS BE N?
  "      DcolLocal < N;\n"
  "      Dcol += get_global_size(1),DcolLocal += get_global_size(1)) {\n"

      // Initialise the accumulation register
  "      acc = 0.0;\n"
   if (need_transpose) {
  // inner loop
    result +=
  "      for(Dinner=0,DinnerA = DmatrixA + DrowLocal,DinnerB = DmatrixB + DcolLocal;"
  "        Dinner < K;\n"
  "        Dinner++,DinnerA += NpadA,DinnerB += NpadB);\n"
  //  "        acc+= A[DrowA + Dinner * NpadA] * B[DcolB + Dinner * NpadB];\n "     
  // do cache
  // entry C[Drow, Dcol] = inner(A[, Drow], B[, Dcol])
  // copy A[seq(DrowLocal, len= get_local_size(0) ), Drow] to Acache
  // DrowLocal is Dinner
  "wait = async_work_group_copy("
  "  Acache,"
  "  A[DinnerA]," //
  "  get_local_size(0),"
  "  0);\n"

     } else{ // not transposed
 // inner loop
  "      for(Dinner=0,DinnerA = DmatrixA + DrowLocal*NpadA, DinnerB = DmatrixB + DcolLocal;"
  "        Dinner < K;\n"
  "        Dinner++, DinnerA++, DinnerB += NpadB);\n"
  //  "        acc+= A[DrowA + Dinner * NpadA] * B[DcolB + Dinner * NpadB];\n "     
  // do cache
  // copy A[,Drow]
  "wait = async_work_group_strided_copy("
  "  Acache,"
  "  A[DinnerA],"
  "  get_local_size(0),"
  "  NpadA,"
  "  0);\n"
}
// copy B[Dcol, ] to Bcache
  result +=
  "wait = async_work_group_copy("
  "  Bcache,"
  "  B[DinnerB],"
  "  get_local_size(1),"
  "  wait);\n";
  "wait_group_events(1, &wait);\n"

  "        acc+= Acache[get_local_id(0)] * Bcache[get_local_id(1)];\n "     
  "      }//Dinner\n";

result +=  "C[DrowC0 + Dcol] = acc;\n"
     "}\n" //Dcol
     "}\n"  	//Drow
     "}\n"//Dmatrix
     
     "}\n"; // function

  return(result);
}













template <typename T> 
void gemmBatch(
    viennacl::matrix<T> &A,
    viennacl::matrix<T> &B,
    viennacl::matrix<T> &C,
    const int z, //num of batches in row//std::vector<int> Nglobal,//   std::vector<int> Nlocal,// const int NlocalCache, 
    const int zc, //num of batches in col
    const int need_transpose,
    const IntegerVector Nglobal,
    const int ctx_id) {
  
  int M,K;
  const int N = B.size2()/zc;
  std::string gemmString;
 
  viennacl::ocl::context ctx(viennacl::ocl::get_context(ctx_id));
  
  cl_device_type type_check = ctx.current_device().type();
  

  if (need_transpose==1){
    M = A.size2()/zc;
    K = A.size1()/z;
    gemmString =  gemmBatchString<T>(  
      M, 
      N, // ncol of Bi
      K,
      A.internal_size2(), 
      B.internal_size2(), 
      C.internal_size2(),
      z, // 
      B.size2(),
      A.internal_size2()*K,//NpadBetweenMatricesA,
      B.internal_size2()*K,//NpadBetweenMatricesB,
      C.internal_size2()*M,
      need_transpose);
  }
  else if (need_transpose==0){
    M = A.size1()/z;
    K = A.size2()/zc;
    gemmString =  gemmBatchString<T>(  
      M, 
      N, // ncol of Bi
      K,
      A.internal_size2(), 
      B.internal_size2(), 
      C.internal_size2(),
      z, // 
      B.size2(),
      A.internal_size2()*M,//NpadBetweenMatricesA,
      B.internal_size2()*K,//NpadBetweenMatricesB,
      C.internal_size2()*M,
      need_transpose);
  }

  
#ifdef DEBUG
  Rcpp::Rcout << gemmString << "\n\n";
#endif  
  
  
  viennacl::ocl::program & my_prog = ctx.add_program(gemmString, "mymykernel");

  viennacl::ocl::kernel & gemmKernel = my_prog.get_kernel("gemm2");

  gemmKernel.global_work_size(0, Nglobal[0]);
  gemmKernel.global_work_size(1, Nglobal[1]);
  gemmKernel.global_work_size(2, Nglobal[2]);
  
  //gemmKernel.local_work_size(0, Nlocal[0]);
  //gemmKernel.local_work_size(1, Nlocal[1]);
 viennacl::ocl::enqueue(gemmKernel( A, B, C));
}


template <typename T> 
SEXP gemmBatchTyped( Rcpp::S4 AR,
                     Rcpp::S4 BR,
                     Rcpp::S4 CR,
                     const int z,
                     const int zc, //num of batches in col
                     const int need_transpose,
                     Rcpp::IntegerVector NglobalR) {
  
  
  const int ctx_id = INTEGER(CR.slot(".context_index"))[0]-1;
  const bool BisVCL=1;
  
  std::shared_ptr<viennacl::matrix<T> > A = getVCLptr<T>(AR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > B = getVCLptr<T>(BR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > C = getVCLptr<T>(CR.slot("address"), BisVCL, ctx_id);
  
  gemmBatch<T>(*A, *B, *C, z, zc, need_transpose, NglobalR, ctx_id);	
  
  return Rcpp::wrap(0L);
}




//' 
//' Multiplies a rectangular matrix by a rectangular matrix in batches
//'
//' @param C output matrices, stacked row-wise
//' @param A rectangular matrices
//' @param B rectangular matrices 
//' @param Nglobal vector of number of global work items//' @param Nlocal vector of number of local work items//' @param NlocalCache elements in local cache
//' @export
// [[Rcpp::export]]
SEXP gemmBatchBackend(
    Rcpp::S4 A,
    Rcpp::S4 B,
    Rcpp::S4 C,  //output matrices, stacked row-wise
    const int z, //num of batches in row
    const int zc, //num of batches in col
    const int need_transpose, //A need transpose?
    Rcpp::IntegerVector Nglobal) {
  
//#ifdef UNDEF  
  Rcpp::traits::input_parameter< std::string >::type classVarR(RCPP_GET_CLASS(C));
  std::string precision_type = (std::string) classVarR;
  if(precision_type == "fvclMatrix") {
    gemmBatchTyped<float>(A, B, C, z, zc, need_transpose, Nglobal);
  } else if (precision_type == "dvclMatrix") {
    gemmBatchTyped<double>(A, B, C, z, zc, need_transpose, Nglobal);
  } else {
    Rcpp::warning("class of var must be fvclMatrix or dvclMatrix");
  }
//#endif  
    return Rcpp::wrap(0L);

}


****/

















