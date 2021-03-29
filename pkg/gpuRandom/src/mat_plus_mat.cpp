#include "gpuRandom.hpp"


//#define DEBUGKERNEL
template <typename T> 
std::string matrix_plus_matrixString(const int Nrow, const int Ncol, const int NpadCola, const int NpadColb, const int NpadColc) {  //internal column size
  
  std::string typeString = openclTypeString<T>();  // type of the sum of log factorial
  
  std::string result = "";
  
  if(typeString == "double") {
    result += "\n#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n";
  }
  
  
  result +=
    "\n#define Nrow " + std::to_string(Nrow) + "\n"    
    "#define Ncol " + std::to_string(Ncol) + "\n"
    "#define NpadCola " + std::to_string(NpadCola) + "\n"
    "#define NpadColb " + std::to_string(NpadColb) + "\n"
    "#define NpadColc " + std::to_string(NpadColc) + "\n";    
  
  
  result += 
    "\n\n__kernel void matrix_add_matrix(\n"
    "  __global " + typeString + "* a,\n"  
    "  __global " + typeString + "* b,\n"
    "  __global " + typeString + "* result"  
    "){\n\n";  
  
  
  result += "int Drow, Dcol;\n";
  
  
  
    result += 
      "  for(Drow = get_global_id(0);   Drow < Nrow;  Drow += get_global_size(0)){\n"
      "    for(Dcol = get_global_id(1);  Dcol < Ncol;   Dcol += get_global_size(1)){\n"
      
      "  result[Drow*NpadColc+Dcol] = a[Drow*NpadCola+Dcol] + b[Drow*NpadColb] ;\n"
      
      "    } // end loop through columns\n"
      "  } // end loop through rows\n";
    
  
  
  result += 
    "}\n";
  
  
  return(result);
}










template<typename T> 
void matrix_matrix_sum(
    viennacl::matrix<T> &a,//
    viennacl::matrix<T> &b,  // must be a matrix of 1 columnn
    viennacl::matrix<T> &sum,
    Rcpp::IntegerVector numWorkItems,
    int ctx_id) {
  
  if (b.size1() > 1  && (b.size2() > 1) ){
    Rcpp::Rcout << "Error: cannot do plus operation" << "\n\n";
    // return EXIT_FAILURE;
  }
  
  std::string KernelString = matrix_plus_matrixString<T>(
    a.size1(), 
    a.size2(),
    a.internal_size2(),
    b.internal_size2(),
    sum.internal_size2()
  );
  
  viennacl::ocl::switch_context(ctx_id);
  viennacl::ocl::program & my_prog = viennacl::ocl::current_context().add_program(KernelString, "my_kernel");
  
#ifdef DEBUGKERNEL
  Rcpp::Rcout << KernelString << "\n\n";
#endif  
  
  viennacl::ocl::kernel &matrix_plus_matrixKernel = my_prog.get_kernel("matrix_add_matrix");
  matrix_plus_matrixKernel.global_work_size(0, numWorkItems[0]);
  matrix_plus_matrixKernel.global_work_size(1, numWorkItems[1]);
  
  viennacl::ocl::enqueue(matrix_plus_matrixKernel(a, b, sum) );
  
  
  // return 1L;
  
}



template<typename T> 
void matrix_matrix_sumTemplated(
    Rcpp::S4 aR,
    Rcpp::S4 bR,
    Rcpp::S4 sumR,
    Rcpp::IntegerVector numWorkItems) {
  
  
  const bool BisVCL=1;
  const int ctx_id = INTEGER(sumR.slot(".context_index"))[0]-1;
  std::shared_ptr<viennacl::matrix<T> > a = getVCLptr<T>(aR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > b = getVCLptr<T>(bR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > sum = getVCLptr<T>(sumR.slot("address"), BisVCL, ctx_id);
  
  matrix_matrix_sum(*a, *b, *sum, numWorkItems, ctx_id);
  
}





//[[Rcpp::export]]
void matrix_matrix_sumBackend(
    Rcpp::S4 aR,
    Rcpp::S4 bR,
    Rcpp::S4 sumR,
    Rcpp::IntegerVector numWorkItems) {
  
  
  Rcpp::traits::input_parameter< std::string >::type classVarR(RCPP_GET_CLASS(sumR));
  std::string precision_type = (std::string) classVarR;
  
  
  if(precision_type == "fvclMatrix") {
    return (matrix_matrix_sumTemplated<float>(aR, bR, sumR, numWorkItems));
  } else if (precision_type == "dvclMatrix") {
    return (matrix_matrix_sumTemplated<double>(aR, bR, sumR, numWorkItems));
  } else if (precision_type  == "ivclMatrix") {
    return( matrix_matrix_sumTemplated<int>(aR, bR, sumR, numWorkItems));
  }
  
}








