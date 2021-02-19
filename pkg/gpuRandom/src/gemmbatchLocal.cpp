

/*
Rcpp::IntegerVector workgroupSize,   // global 0 1 2, local 0 1 2
Rcpp::IntegerVector transposeABC, // transposeA, transposeB, transposeC
Rcpp::IntegerVector submatrixA, // [rowStart, nRowsSub, nRowsTotal, colStart, ...]
Rcpp::IntegerVector submatrixB,
Rcpp::IntegerVector submatrixC,
Rcpp::IntegerVector batches, // nRow, nCol, recycleArow, recycleAcol, recycleB row col
Rcpp::IntegerVector NlocalCache,  //cacheSizeA, cacheSizeB, 
 
NOTE: recycling A, B not implemented 
 */

/*
  * submatrixA = 
  *  [rowStartA, nRowsAsub, nRowsA, colStartA, nColsAsub, nColsA]
  *  indexing from zero
  *  submatrix[DmatrixRow, DmatrixCol] of A is
  *  A[seq(nRowsA * DmatrixRow + rowStartA, len=nRowsAsub),
  *  seq(nColsA * DmatrixCol + colStartA, len=nColsAsub)]
  *  
  * batches = 
  * nRowBatch, nColBatch, recycleArow, recycleAcol, recycleBrow, recycleBcol
  * recycleArow = 1, there are no row batches for A, use the same A for all batches
  * 
  * Rcpp::IntegerVector workgroupSize, global 0 1 2, local 0 1 2
  *  matrix global0 rowbatch local0 colbatch, rows of C (and A), cols of C (and B)
  * DmatrixRow = global0; Dcolbatch = local0, Drow = global1, Dcol = global2
  * 
  * 
  * A[DmatrixRow, DmatrixCol][Drow, Dcol] =
  * A[DmatrixRow * (NpadA * nRowsTotal) + 
  *     NpadA * (rowStart + Drow) + 
  *     DmatrixCol * nColsTotal +
  *     Dcol]
  * 
  */
 
#include "gpuRandom.hpp"
using namespace Rcpp;


//#define DEBUG


template <typename T> 
std::string gemmBatch2String(
    Rcpp::IntegerVector transposeABC,  
    Rcpp::IntegerVector submatrixA,  
    Rcpp::IntegerVector submatrixB,
    Rcpp::IntegerVector submatrixC,
    Rcpp::IntegerVector batches,  
    Rcpp::IntegerVector NlocalCache,
    int NpadA, int NpadB, int NpadC) { 

  std::string typeString = openclTypeString<T>();
  std::string result = "";
  

  if(typeString == "double") {
    result += "\n#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n\n";
  }
  
  result += 
    "#define NpadA "+ std::to_string(NpadA) + "\n"  
    "#define NpadB "+ std::to_string(NpadB) + "\n"  
    "#define NpadC "+ std::to_string(NpadC) + "\n\n";  
    
  if(transposeABC[2]) { // transpose C, Nrow is rows of C^T, cols of C
    result += 
      "#define Nrow " + std::to_string(submatrixC[4]) + "\n"  
      "#define rowStartC " + std::to_string(submatrixC[3]) + "\n"
      "#define NrowTotalC " + std::to_string(submatrixC[5]) + "\n"
      "#define Ncol " + std::to_string(submatrixC[1]) + "\n"
      "#define colStartC " + std::to_string(submatrixC[0]) + "\n"
      "#define NcolTotalC " + std::to_string(submatrixC[2]) + "\n"
      "#define Cij = ( DmatrixRow * NpadC * NcolTotalC + NpadC * (colStartC + Dcol) + DmatrixCol * NrowTotalC + rowStartC + Drow)\n\n";
  } else {
      result += 
        "#define Nrow " + std::to_string(submatrixC[1]) + "\n"  
        "#define rowStartC " + std::to_string(submatrixC[0]) + "\n"
        "#define NrowTotalC " + std::to_string(submatrixC[2]) + "\n"
        "#define Ncol " + std::to_string(submatrixC[4]) + "\n"
        "#define colStartC " + std::to_string(submatrixC[3]) + "\n"
        "#define NcolTotalC " + std::to_string(submatrixC[5]) + "\n"
        "#define Cij ( DmatrixRow * NpadC * NrowTotalC + NpadC * (rowStartC + Drow) + DmatrixCol * NcolTotalC + colStartC + Dcol)\n\n";
  }
  if(transposeABC[0]) { // transpose A, Ninner is cols of A^T, rows of A
    result += 
    "#define Ninner " + std::to_string(submatrixA[4]) + "\n"
    "#define rowStartA " + std::to_string(submatrixA[3]) + "\n"
    "#define NrowTotalA " + std::to_string(submatrixA[5]) + "\n"
    "#define colStartA " + std::to_string(submatrixA[0]) + "\n"
    "#define NcolTotalA " + std::to_string(submatrixA[2]) + "\n"
    "#define Ai0(ii) ( DmatrixRow * NpadA * NcolTotalA + NpadA * (colStartA + ii) + DmatrixCol * NrowTotalA + rowStartA)\n\n";
  } else { // no transpose, Ninner is cols of A
    result += 
      "#define Ninner " + std::to_string(submatrixA[1]) + "\n"
    "#define rowStartA " + std::to_string(submatrixA[0]) + "\n"
    "#define NrowTotalA " + std::to_string(submatrixA[2]) + "\n"
    "#define colStartA " + std::to_string(submatrixA[3]) + "\n"
    "#define NcolTotalA " + std::to_string(submatrixA[5]) + "\n"
    "#define Ai0(ii) ( DmatrixRow * NpadA * NrowTotalA + NpadA * (rowStartA + ii) + DmatrixCol * NcolTotalA + colStartA)\n\n";
  }

  if(transposeABC[1]) { // transpose B
    result += 
      "#define rowStartB " + std::to_string(submatrixA[3]) + "\n"
      "#define NrowTotalB " + std::to_string(submatrixA[5]) + "\n"
      "#define colStartB " + std::to_string(submatrixA[0]) + "\n"
      "#define NcolTotalB " + std::to_string(submatrixA[2]) + "\n"
    "#define Bi0(ii) ( DmatrixRow * NpadB * NcolTotalB + NpadB * (colStartB + ii) + DmatrixCol * NrowTotalB + rowStartB)\n\n";
  } else { // no transpose, Ninner is cols of A
    result += 
    "#define rowStartB " + std::to_string(submatrixB[0]) + "\n"
    "#define NrowTotalB " + std::to_string(submatrixB[2]) + "\n"
    "#define colStartB " + std::to_string(submatrixB[3]) + "\n"
    "#define NcolTotalB " + std::to_string(submatrixB[5]) + "\n"
    "#define Bi0(ii) ( DmatrixRow * NpadB * NrowTotalB + NpadB * (rowStartB + ii) + DmatrixCol * NcolTotalB + colStartB)\n\n";
  }

    result +=  
    "#define NmatrixRow " + std::to_string(batches[0]) + "\n"
    "#define NmatrixCol " + std::to_string(batches[1]) + "\n"; 
  
  result +=  
  "#define cacheSizeA " + std::to_string(NlocalCache[0]) + "\n"
  "#define cacheSizeB " + std::to_string(NlocalCache[1]) + "\n"; 
  
    

  result += "\n__kernel void gemm( __global "  + typeString+ "* A,\n"
                               " __global "  + typeString+ "* B,\n"
                               " __global " + typeString+ "* C){\n\n";
                              
  
  result += typeString + " acc;\n"
           "int DmatrixRow, DmatrixCol;\n"
           "int DthisA, DthisB, DthisC;\n"
           "int Drow, DrowBlock, Dcol, DcolBlock;\n"
           "int Dinner, Dinnerp1, Dinnerp2;\n";
  result += "event_t wait;\n\n";

  result += "local " + typeString + 
    " localCacheA1[cacheSizeA], localCacheB1[cacheSizeB];\n";   
    result += "local " + typeString + 
      " localCacheA2[cacheSizeA], localCacheB2[cacheSizeB];\n";

result += "\n";
result += 
  " for(DmatrixRow = get_group_id(0);\n" 
  "   DmatrixRow < NmatrixRow;\n" 
  "   DmatrixRow += get_num_groups(0)) {\n"
  " for(DmatrixCol = get_local_id(1);\n"
  "   DmatrixCol < NmatrixCol;\n" 
  "   DmatrixCol += get_local_size(1)) {\n\n";
result +=
  "   for(DrowBlock = get_group_id(1)*get_local_size(1);\n" // row for work item ?,0,?
  "     DrowBlock < Nrow;\n"
  "     DrowBlock += get_global_size(1)) {\n"
  "   Drow = DrowBlock + get_local_id(1);\n";

  result += "\n/* cache A keep */\n\n";


result +=
  "   for(DcolBlock = get_group_id(2)*get_local_size(2);\n" // col for work item ?,?,0
  "     DcolBlock < Ncol;\n"
  "     DcolBlock += get_global_size(2)) {\n"
  "   Dcol = DcolBlock + get_local_id(2);\n";
  
result += "acc = 0.0;\n";

result += "wait =  (event_t) 0;\n";

if(transposeABC[0]) { 
  result +=  "  wait = async_work_group_copy(\n"
  "    localCacheA1, &A[Ai0(0)],\n"
  "    get_local_size(1), (event_t) 0);\n";
} else {
  result +=  "  wait = async_work_group_strided_copy(\n"
  "    localCacheA1, &A[Ai0(0)],\n"
  "    get_local_size(1), NpadA, (event_t) 0);\n";
}
if(transposeABC[1]) { 
  result +=  "  wait = async_work_group_strided_copy(\n"
  "    localCacheB1, &B[Bi0(0)],\n"
  "    get_local_size(2), NpadB, wait);\n";
} else {
  result +=  "  wait = async_work_group_copy(\n"
  "    localCacheB1, &B[Bi0(0)],\n"
  "    get_local_size(2), wait);\n";
}

  result += // loop through remaining inner
    "   for(Dinner = 0,Dinnerp1=1,Dinnerp2=2;\n" 
    "     Dinner < Ninner;\n"
    "     Dinner++, Dinnerp1++, Dinnerp2++) {\n";
  // cache row Dinner of C, col Dinner of A
  if(transposeABC[0]) { 
    result +=  "  wait = async_work_group_copy(\n"
    "    localCacheA1, &A[Ai0(Dinner)],\n"
    "    get_local_size(1), (event_t) 0);\n";
  } else {
result +=  "  wait = async_work_group_strided_copy(\n"
  "    localCacheA1, &A[Ai0(Dinner)],\n"
  "    get_local_size(1), NpadA, (event_t) 0);\n";
  }
  if(transposeABC[1]) { 
    result +=  "  wait = async_work_group_strided_copy(\n"
    "    localCacheB1, &B[Bi0(Dinner)],\n"
    "    get_local_size(2), NpadB, wait);\n";
  } else {
    result +=  "  wait = async_work_group_copy(\n"
    "    localCacheB1, &B[Bi0(Dinner)],\n"
    "    get_local_size(2), wait);\n";
  }
  result += "  wait_group_events (1, &wait);\n";
//  "  barrier(CLK_LOCAL_MEM_FENCE);\n";

  result += "\n";
  result += "acc += localCacheA1[get_local_id(1)] * localCacheB1[get_local_id(2)];\n";
  result += "\n";

  result += 
    "   }//Dinner Ninner\n";
  
  result += 
    "\nC[Cij] = acc;\n\n";
    

result += 
  "   }//DrowBlock\n"
  "   }//DcolBlock\n";
   
   result += 
  " }//DmatrixRow\n"
  " }//DmatrixCol\n";
   result += 
     " }//kernel\n";
     
  
  return(result);
}








template <typename T> 
int gemmBatch2(
    viennacl::matrix<T> &A,
    viennacl::matrix<T> &B,
    viennacl::matrix<T> &C,
    Rcpp::IntegerVector transposeABC,  
    Rcpp::IntegerVector submatrixA,
    Rcpp::IntegerVector submatrixB,
    Rcpp::IntegerVector submatrixC,  
    Rcpp::IntegerVector batches, 
    Rcpp::IntegerVector workgroupSize,
    Rcpp::IntegerVector NlocalCache, 
    const int verbose,
    const int ctx_id) {
  

  std::string gemmString = gemmBatch2String<T>(
    transposeABC,  
    submatrixA,
    submatrixB,
    submatrixC,
    batches, 
    NlocalCache, 
    (int) A.internal_size2(),
    (int) B.internal_size2(),
    (int) C.internal_size2()
  );
  
  if(verbose) Rcpp::Rcout << "\n\n" << gemmString << "\n\n";

  viennacl::ocl::context ctx(viennacl::ocl::get_context(ctx_id));
  viennacl::ocl::program & my_prog = ctx.add_program(gemmString, "mykernel");
  viennacl::ocl::kernel & gemmKernel = my_prog.get_kernel("gemm");
  
  gemmKernel.global_work_size(0, workgroupSize[0]);
  gemmKernel.global_work_size(1, workgroupSize[1]);
  gemmKernel.global_work_size(2, workgroupSize[2]);
  gemmKernel.local_work_size(0, workgroupSize[3]);
  gemmKernel.local_work_size(1, workgroupSize[4]);
  gemmKernel.local_work_size(2, workgroupSize[5]);  
  
  viennacl::ocl::enqueue(gemmKernel(A, B, C));

  return 1;
}


template <typename T> 
SEXP gemmBatch2Typed(Rcpp::S4 AR,
                     Rcpp::S4 BR,
                     Rcpp::S4 CR,
                     Rcpp::IntegerVector transposeABC,  
                     Rcpp::IntegerVector submatrixA,
                     Rcpp::IntegerVector submatrixB,
                     Rcpp::IntegerVector submatrixC,  
                     Rcpp::IntegerVector batches, 
                     Rcpp::IntegerVector workgroupSize,
                     Rcpp::IntegerVector NlocalCache, 
                     const int verbose) {
  
  
  const int ctx_id = INTEGER(CR.slot(".context_index"))[0]-1;
  const bool BisVCL=1;
  int result;
  
  std::shared_ptr<viennacl::matrix<T> > A = getVCLptr<T>(AR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > B = getVCLptr<T>(BR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > C = getVCLptr<T>(CR.slot("address"), BisVCL, ctx_id);
  
  result = gemmBatch2<T>(*A, *B, *C, transposeABC, 
                submatrixA, submatrixB, submatrixC, batches, 
                workgroupSize, NlocalCache, verbose, ctx_id);	
  
  return Rcpp::wrap(result);
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
SEXP gemmBatch2backend(
    Rcpp::S4 A,
    Rcpp::S4 B,  
    Rcpp::S4 C,
    Rcpp::IntegerVector transposeABC,  
    Rcpp::IntegerVector submatrixA,
    Rcpp::IntegerVector submatrixB,
    Rcpp::IntegerVector submatrixC, 
    Rcpp::IntegerVector batches, 
    Rcpp::IntegerVector workgroupSize,   
    Rcpp::IntegerVector NlocalCache,
    const int verbose) {
  
  SEXP result;
  
//#ifdef UNDEF  
  Rcpp::traits::input_parameter< std::string >::type classVarR(RCPP_GET_CLASS(C));
  std::string precision_type = (std::string) classVarR;
  if(precision_type == "fvclMatrix") {
    result=gemmBatch2Typed<float>(A, B, C, 
                                  transposeABC, 
                           submatrixA, submatrixB, submatrixC, batches, 
                           workgroupSize, NlocalCache, verbose);
  } else if (precision_type == "dvclMatrix") {
    result=gemmBatch2Typed<double>(A, B, C, 
                                   transposeABC, 
                            submatrixA, submatrixB, submatrixC, batches, 
                            workgroupSize, NlocalCache, verbose);
  } else {
    Rcpp::warning("class of var must be fvclMatrix or dvclMatrix");
  }
//#endif  
    return result;

}




















