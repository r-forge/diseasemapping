
#include "gpuRandom.hpp"
using namespace Rcpp;


//#define DEBUG

// C = A B   C_i is M by N, A_i is M by K, B_i is K by N  (if no transposition)
// work items w1, w2, w3, i.e 32 by 32 by 4
// for now local items 1x1x1
// suppose matrices C are 120 by 140, 10 matrices
// item w1, w2, 0  will do matrices 0,4,8.  
// item w1, w2, 3 will do matrices 3, 7

template <typename T> 
std::string gemmBatchString(
    const int M, // rows of Ai
    const int N, // columns of Bi
    const int K, // common size of Ai Bi
    const int rowbatch, //number of batches in row
    const int Acolbatch, 
    const int Bcolbatch,
    const int NpadA, 
    const int NpadB, 
    const int NpadC,
    const int NpadMatrixA, 
    const int NpadMatrixB, 
    const int NpadMatrixC,
    const int need_transpose) { 
  
 
  std::string typeString = openclTypeString<T>();
  std::string result = "";
  
  if(typeString == "double") {
    result += "\n#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n\n";
  }
  
  result += 
    "#define M " + std::to_string(M) + "\n" // rows 
    "#define N " + std::to_string(N) + "\n" // columns
    "#define K " + std::to_string(K) + "\n" // inter
    "#define rowbatch " + std::to_string(rowbatch) + "\n" // matrices
    "#define Acolbatch " + std::to_string(Acolbatch) + "\n" 
    "#define Bcolbatch " + std::to_string(Bcolbatch) + "\n" //num of batches in B columnwise
    "#define NpadA " + std::to_string(NpadA) + "\n"    
    "#define NpadB " + std::to_string(NpadB) + "\n"  
    "#define NpadC " + std::to_string(NpadC) + "\n"
    "#define NpadMatrixA " + std::to_string(NpadMatrixA) + "\n"    
    "#define NpadMatrixB " + std::to_string(NpadMatrixB) + "\n"  
    "#define NpadMatrixC " + std::to_string(NpadMatrixC) + "\n"
    "#define need_transpose  " + std::to_string(need_transpose) + "\n";
  
  result += "__kernel void gemm2( __global "  + typeString+ "* A,\n"
                                " __global "  + typeString+ "* B,\n"
                                " __global " + typeString+ "* C){\n";
                              
      
      
      // Initialise the accumulation register
      //#define Ncache 100	
    
  result += typeString + " acc ;\n"
         "int numbatch, Dmatrix, Drow, DrowA0, DrowA, Dcol, DrowC0, Dinner, DrowB0, DcolB;\n"
         "int Ccoltotal= N * Acolbatch;\n";
      
  result += "for (Dmatrix=get_global_id(2); Dmatrix < rowbatch; Dmatrix += get_global_size(2)) {\n"
                "DrowB0 = Dmatrix * NpadMatrixB;\n"   
                "for (Drow = get_global_id(0); Drow < M; Drow += get_global_size(0)) {\n"
                "DrowC0 = Dmatrix * NpadMatrixC + Drow*NpadC;\n";
          
if (need_transpose && (Bcolbatch == Acolbatch)) {
          
result += "DrowA0 = Dmatrix * NpadMatrixA + Drow;\n"
            
            "for (Dcol = get_global_id(1); Dcol < Ccoltotal; Dcol += get_global_size(1)) {\n"
              "DcolB = DrowB0 + Dcol;\n"
              "acc = 0;\n"
              "numbatch = Dcol/N;\n"
              "DrowA = DrowA0 + M*numbatch;\n"
              "for(Dinner = 0; Dinner < K; Dinner++){\n "
                "acc+= A[DrowA + Dinner * NpadA] * B[DcolB + Dinner * NpadB];\n "     
             " }\n";
                 
}else if (need_transpose && (Bcolbatch == 1) && (Acolbatch > 1)){
  
result +=  " DrowA0 = Dmatrix * NpadMatrixA + Drow;\n" 
       " for (Dcol = get_global_id(1); Dcol < Ccoltotal; Dcol += get_global_size(1)) {\n"
    
    "acc = 0;\n"
    "numbatch=Dcol/N;\n"
    "DrowA = DrowA0 + M*numbatch;\n"
   " DcolB = DrowB0 + Dcol - numbatch*N;\n"
    
   " for(Dinner = 0; Dinner < K; Dinner++){\n "
     " acc+= A[DrowA + Dinner * NpadA] * B[DcolB + Dinner * NpadB];\n  "   
   " }\n";
}else if (!need_transpose && (Bcolbatch == Acolbatch)){
  
result += 
     "DrowA0 = Dmatrix * NpadMatrixA + Drow*NpadA;\n"  // points to entry (Drow,0) of submatrix Dmatrix of A
     
     "for (Dcol = get_global_id(1); Dcol < Ccoltotal; Dcol += get_global_size(1)) {\n"
     "DcolB = DrowB0 + Dcol;\n"            
     "numbatch = Dcol/N;\n"
     "DrowA = DrowA0 + K*numbatch;\n"
     "acc = 0;\n"
     "for(Dinner = 0; Dinner < K; Dinner++){\n"
     "acc+= A[DrowA + Dinner] * B[DcolB + Dinner * NpadB];\n"
     " }\n";
}
else if (!need_transpose && (Bcolbatch == 1) && (Acolbatch > 1)){
result += "DrowA0 = Dmatrix * NpadMatrixA + Drow * NpadA;\n"     
  "for (Dcol = get_global_id(1); Dcol < Ccoltotal; Dcol += get_global_size(1)) {\n "  
    
    "numbatch=Dcol/N;\n"
    "DrowA = DrowA0 + K*numbatch;\n"
    "DcolB = DrowB0 + Dcol-numbatch*N;\n"
    "acc = 0;\n"
   " for(Dinner = 0; Dinner < K; Dinner++){\n "
      "acc+= A[DrowA+ Dinner] * B[DcolB + Dinner * NpadB];\n"
    "}\n";
}

result +=  "C[DrowC0 + Dcol] = acc;\n"
     "}\n" //Dcol
     "}\n"  	//Drow
     "}\n"//Dmatrix
     
     "}\n";
     
  return(result);
}













template <typename T> 
void gemmBatch(
    viennacl::matrix<T> &A,
    viennacl::matrix<T> &B,
    viennacl::matrix<T> &C,
    const int rowbatch, //num of batches in row//std::vector<int> Nglobal,//   std::vector<int> Nlocal,// const int NlocalCache, 
    const int Acolbatch, //num of batches in col
    const int Bcolbatch,
    const int need_transpose,
    const IntegerVector Nglobal,
    const int ctx_id) {
  
  int M,K;
  const int N = B.size2()/Bcolbatch;
  std::string gemmString;
 
  viennacl::ocl::context ctx(viennacl::ocl::get_context(ctx_id));
  
  cl_device_type type_check = ctx.current_device().type();
  
 
  if (need_transpose==1){
    M = A.size2()/Acolbatch;
    K = A.size1()/rowbatch;
    gemmString =  gemmBatchString<T>(  
      M, 
      N, // ncol of Bi
      K,
      rowbatch,
      Acolbatch,
      Bcolbatch,
      A.internal_size2(), 
      B.internal_size2(), 
      C.internal_size2(),
      A.internal_size2()*K,//NpadBetweenMatricesA,
      B.internal_size2()*K,//NpadBetweenMatricesB,
      C.internal_size2()*M,
      need_transpose);
  }
  else if (need_transpose==0){
    M = A.size1()/rowbatch;
    K = A.size2()/Acolbatch;
    gemmString =  gemmBatchString<T>(  
      M, 
      N, // ncol of Bi
      K,
      rowbatch,
      Acolbatch,
      Bcolbatch,
      A.internal_size2(), 
      B.internal_size2(), 
      C.internal_size2(),
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
                     const int rowbatch, //num of batches in row
                     const int Acolbatch, //num of batches in col
                     const int Bcolbatch,
                     const int need_transpose,
                     Rcpp::IntegerVector NglobalR) {
  
  
  const int ctx_id = INTEGER(CR.slot(".context_index"))[0]-1;
  const bool BisVCL=1;
  
  std::shared_ptr<viennacl::matrix<T> > A = getVCLptr<T>(AR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > B = getVCLptr<T>(BR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > C = getVCLptr<T>(CR.slot("address"), BisVCL, ctx_id);
  
  gemmBatch<T>(*A, *B, *C, rowbatch, Acolbatch, Bcolbatch, need_transpose, NglobalR, ctx_id);	
  
  return Rcpp::wrap(0L);
}





// [[Rcpp::export]]
SEXP gemmBatchBackend(
    Rcpp::S4 A,
    Rcpp::S4 B,
    Rcpp::S4 C,  //output matrices, stacked row-wise
    const int rowbatch, //num of batches in row
    const int Acolbatch, //num of batches in col
    const int Bcolbatch,
    const int need_transpose, //A need transpose?
    Rcpp::IntegerVector Nglobal) {
  
//#ifdef UNDEF  
  Rcpp::traits::input_parameter< std::string >::type classVarR(RCPP_GET_CLASS(C));
  std::string precision_type = (std::string) classVarR;
  if(precision_type == "fvclMatrix") {
    gemmBatchTyped<float>(A, B, C, rowbatch, Acolbatch, Bcolbatch, need_transpose, Nglobal);
  } else if (precision_type == "dvclMatrix") {
    gemmBatchTyped<double>(A, B, C, rowbatch, Acolbatch, Bcolbatch, need_transpose, Nglobal);
  } else {
    Rcpp::warning("class of var must be fvclMatrix or dvclMatrix");
  }
//#endif  
    return Rcpp::wrap(0L);

}




















