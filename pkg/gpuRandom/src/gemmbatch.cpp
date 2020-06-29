
#include "gpuRandom.hpp"
using namespace Rcpp;


#define DEBUG

// C = A B   C_i is M by N, A_i is M by K B_i is K by N 
// work items w1, w2, w3 global size, i.e 32 by 32 by 4
// for now local items 1x1x1
// suppose matrices C are 120 by 140, 10 matrices
// item w1, w2, 0  will do matrices 0,4,8.  
// item w1, w2, 3 will do matrices 3, 7

template <typename T> 
std::string gemmBatchString(
    const int M, // rows of Ai
    const int N, // columns of Bi
    const int K, // common size of Ai Bi
    const int NpadA, 
    const int NpadB, 
    const int NpadC,
    const int z, //number of batches
    const int NpadMatrixA, 
    const int NpadMatrixB, 
    const int NpadMatrixC) { 
  
 
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
    "#define NpadA " + std::to_string(NpadA) + "\n"    
    "#define NpadB " + std::to_string(NpadB) + "\n"  
    "#define NpadC " + std::to_string(NpadC) + "\n"
    "#define NpadMatrixA " + std::to_string(NpadMatrixA) + "\n"    
    "#define NpadMatrixB " + std::to_string(NpadMatrixB) + "\n"  
    "#define NpadMatrixC " + std::to_string(NpadMatrixC) + "\n";
  
  result += "__kernel void gemm2( __global"  + typeString+ "* A,\n"
                               " __global"  + typeString+ "* B,\n"
                               " __global" + typeString+ "* C) {\n";
      
      
      // Initialise the accumulation register
      //#define Ncache 100	
    
  result += typeString + " acc ;\n"
         "int Dmatrix, Drow, DrowA0, Dcol, DrowC0, Dinner, DrowB, DcolB;\n";
      
      //local float Acache[Ncachje];
#ifdef UNDEF      

      // Loop over all batches
      "for (Dmatrix=get_global_id(2); Dmatrix < z; Dmatrix += get_global_size(2)) {\n"

        " DrowB = Dmatrix * NpadMatrixB;\n"      // B_Dmatrix(0,0)

       " for (Drow = get_global_id(0); Drow < M; Drow += get_global_size(0)) {\n"

          "DrowA0 = Dmatrix * NpadMatrixA + Drow*NpadA;\n" 
          // points to entry (Drow,0) of submatrix Dmatrix of A
          "DrowC0 = Dmatrix * NpadMatrixC + Drow*NpadC;\n" 

          "for (Dcol = get_global_id(1); Dcol < N; Dcol += get_global_size(1)) {\n"
            "DcolB = DrowB + Dcol;\n"            
            "acc = 0;\n"
            
            "for(Dinner = 0; Dinner < K; Dinner++){\n"
              "acc+= A[DrowA0 + Dinner] * B[DcolB + Dinner * NpadB];\n"
           " }\n"
            
            "C[DrowC0 + Dcol] = acc;\n"
          "}\n" //Dcol
        "}\n"  	//Drow
      "}\n"；  //Dmatrix
#endif     
   result +=  "}\n";

  return(result);
}













template <typename T> 
void gemmBatch(
    viennacl::matrix<T> &A,
    viennacl::matrix<T> &B,
    viennacl::matrix<T> &C,
    const int z, //num of batches//std::vector<int> Nglobal,//   std::vector<int> Nlocal,// const int NlocalCache, 
    const IntegerVector Nglobal,
    const int ctx_id) {
  
  const int M = A.size1()/z;
  const int N = B.size2();
  const int K = A.size2();

  // the context
  viennacl::ocl::context ctx(viennacl::ocl::get_context(ctx_id));
  
  cl_device_type type_check = ctx.current_device().type();
  
  std::string gemmString =  gemmBatchString<T>(  
    M, 
    N, // ncol of B
    K,
    A.internal_size2(), 
    B.internal_size2(), 
    C.internal_size2(),
    z, // 
    A.internal_size2()*M,//NpadBetweenMatricesA,
    B.internal_size2()*K,//NpadBetweenMatricesB,
    C.internal_size2()*M);
  
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
  #ifdef UNDEF   
 viennacl::ocl::enqueue(gemmKernel( A, B, C));
#endif 
}


template <typename T> 
SEXP gemmBatchTyped( Rcpp::S4 AR,
                     Rcpp::S4 BR,
                     Rcpp::S4 CR,
                     const int z,
                     Rcpp::IntegerVector NglobalR) {
  
  //std::vector<int> Nglobal = Rcpp::as<std::vector<int> >(NglobalR);
  //std::vector<int> Nlocal = Rcpp::as<std::vector<int> >(NlocalR);
  
  const int ctx_id = INTEGER(CR.slot(".context_index"))[0]-1;
  const bool BisVCL=1;
  
  std::shared_ptr<viennacl::matrix<T> > A = getVCLptr<T>(AR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > B = getVCLptr<T>(BR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > C = getVCLptr<T>(CR.slot("address"), BisVCL, ctx_id);
  
  gemmBatch<T>(*A, *B, *C, z, NglobalR, ctx_id);	
  
  return Rcpp::wrap(0L);
}




//' 
//' Multiplies a rectangular matrix by a rectangular matrix
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
    Rcpp::S4 C,
    const int z,
    Rcpp::IntegerVector Nglobal) {
  
//#ifdef UNDEF  
  Rcpp::traits::input_parameter< std::string >::type classVarR(RCPP_GET_CLASS(C));
  std::string precision_type = (std::string) classVarR;
  if(precision_type == "fvclMatrix") {
    gemmBatchTyped<float>(A, B, C, z, Nglobal);
  } else if (precision_type == "dvclMatrix") {
    gemmBatchTyped<double>(A, B, C, z, Nglobal);
  } else {
    Rcpp::warning("class of var must be fvclMatrix or dvclMatrix");
  }
//#endif  
    return Rcpp::wrap(0L);

}




















