
#include "gpuRandom.hpp"


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
    "#define M " + std::to_string(M) + "\n"    
    "#define N " + std::to_string(N) + "\n"    
    "#define K " + std::to_string(K) + "\n"    
    "#define z " + std::to_string(z) + "\n"  
    "#define NpadA " + std::to_string(NpadA) + "\n"    
    "#define NpadB " + std::to_string(NpadB) + "\n"  
    "#define NpadC " + std::to_string(NpadC) + "\n";
  
  result +=
    "__kernel void GEMM(const __global"  + typeString+ "* A,\n"
                       "const __global"  + typeString+ "* B,\n"
                       "__global" + typeString+ "* C) {\n";
      
      // Thread identifiers
      //  const int row = get_local_id(0); // Local row ID 
      //  const int col = get_local_id(1); // Local col ID 
      
      // Identify the row and column of the C matrix to work on
      // const int globalRow = M*get_group_id(0) + localrow; // Row ID of C 
      // const int globalCol = localcol; // Col ID of C 
      
      // Local memory to fit a tile of TS*TS elements of A and B
      // __local float Asub[M][K];
      // __local float Bsub[K][N];
      
      // Initialise the accumulation register
      //#define Ncache 100	
  result += typeString + " acc = 0.0f;\n"
         "int Dmatrix, Drow, Dcol, Dinner;\n"
      
      //local float Acache[Ncachje];
       "int rowtotal=M*z;\n"
      
      // Loop over all batches
      //const int numBatches = z;
      "for (Dmatrix=get_global_id(2); Dmatrix < z; Dmatrix += get_global_size(2)) {\n"
        
       " for (Drow = Dmatrix*M + get_glboal_id(0); Drow < rowtotal; Drow += get_global_size(0)) {\n"
          //		for(Dcol = 0; Dcol < Ncache; Dcache++) {
          //			Acache[Dcol] = A[Drow * NpadA + Dcol];
          //		}
          "for (Dcol = get_glboal_id(1); Dcol < N; Dcol += get_global_size(1)) {\n"
            
            "acc = 0.0;\n"
            // break stuff
            //for(Dinner = 0; Dinner < Ncache; Dinner++){
            //acc+= Acache[Dinner] * B[Dinner * NpadB + Dcol];
            //}
            // break stuff
            //			for(; Dinner < K; Dinner++){
            "for(Dinner = 0; Dinner < K; Dinner++){\n"
              "acc+= A[Drow*NpadA + Dinner] * B[Dinner * NpadB + Dcol];\n"
           " }\n"
            
            "C[NpadC * Drow + Dcol] = acc;\n"
          "}\n" //Dcol
        "}\n"  	//Drow
      "}\n"  //Dbatch
    "}\n";
  
  return(result);
}













template <typename T> 
void gemmBatch(
    viennacl::matrix<T> &A,
    viennacl::matrix<T> &B,
    viennacl::matrix<T> &C,
    const int z, //num of batches
    std::vector<int> Nglobal,//   std::vector<int> Nlocal,// const int NlocalCache, 
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
    C.internal_size2()*M//NpadBetweenMatricesC,
    );
  
#ifdef DEBUG
  Rcpp::Rcout << clString << "\n\n";
#endif  
  
  
  viennacl::ocl::program & my_prog = ctx.add_program(gemmString, "my_kernel");
  
  viennacl::ocl::kernel & gemmKernel = my_prog.get_kernel("GEMM");
  
  gemmKernel.global_work_size(0, Nglobal[0]);
  gemmKernel.global_work_size(1, Nglobal[1]);
  gemmKernel.global_work_size(2, Nglobal[2]);
  
  //gemmKernel.local_work_size(0, Nlocal[0]);
  //gemmKernel.local_work_size(1, Nlocal[1]);
  
  // diagonals and diagTimesRowOfA
  viennacl::ocl::enqueue(gemmKernel( A, B,C));
  
}





















