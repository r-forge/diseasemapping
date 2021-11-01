#include "gpuRandom.hpp"

//#define DEBUG
//#define DEBUGKERNEL


using namespace Rcpp;
using namespace viennacl; 
using namespace viennacl::linalg;

/*
//#define DEBUGKERNEL
template <typename T> 
std::string LogString(const int Nrow, const int Ncol, const int NpadCol) {  //internal column size
  
  std::string typeString = openclTypeString<T>();

  
  std::string result = "";
  
  if(typeString == "double") {
    result += "\n#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n";
  }
  
  
  result +=
    "\n#define Nrow " + std::to_string(Nrow) + "\n"    
    "#define Ncol " + std::to_string(Ncol) + "\n"
    "#define NpadCol " + std::to_string(NpadCol) + "\n"; 
 
   result += 
   "\n\n__kernel void Logmatrix(\n"
   "  __global " + typeString + "* x,\n"  
   "  __global " + typeString + "* output\n"  
   "){\n\n";  
   
   result += "int Drow, Dcol, Dindex;\n";

     result += 
       
     "  for(Drow = get_global_id(0);   Drow < Nrow;    Drow+=get_global_size(0)){\n"
     "    for(Dcol = get_global_id(1),   Dindex = Drow*NpadCol+Dcol;\n" 
     "        Dcol < Ncol; Dcol+=get_global_size(1), Dindex++){\n"
     "        output[Dindex] =log(x[Dindex]);\n"
     "    } // end loop through columns\n"
     "  } // end loop through rows\n"
     
     
    "}\n";
     
     

  return(result);
}

*/




/*
template <typename T>
void rowsum(
    viennacl::matrix<T> &x,
    viennacl::vector_base<T> &Sum,
    int log,
    Rcpp::IntegerVector numWorkItems,
    int ctx_id) {

    const int nrow = x.size1(), ncol = x.size2(), nPadCol= x.internal_size2();
    
    
    if (log ==1 ){
      viennacl::matrix<T> output(nrow,ncol);  
    
      std::string KernelString = LogString<T>(nrow, ncol, nPadCol);
    // the context
      viennacl::ocl::switch_context(ctx_id);
      viennacl::ocl::program & my_prog = viennacl::ocl::current_context().add_program(KernelString, "my_kernel");
      
#ifdef DEBUGKERNEL
      Rcpp::Rcout << KernelString << "\n\n";
#endif  
      
      
      viennacl::ocl::kernel &LogKernel = my_prog.get_kernel("Logmatrix");
      
      LogKernel.global_work_size(0, numWorkItems[0]);
      LogKernel.global_work_size(1, numWorkItems[1]);
      
      
      viennacl::ocl::enqueue(LogKernel(x, output) );
      
      row_sum_impl(output, Sum);
      
    }else if(log==0){
      
    row_sum_impl(x, Sum);}
}

*/






template <typename T>
void rowsum(
    viennacl::matrix<T> &x,
    viennacl::vector_base<T> &Sum,
    std::string type,
    int log){// Rcpp::IntegerVector numWorkItems,int ctx_id) {
  
  const int nrow = x.size1(), ncol = x.size2();
  
  if (log ==1 ){
    viennacl::matrix<T> logx(nrow,ncol);  
    logx = viennacl::linalg::element_log(x);
    
    if(type=="row"){
    row_sum_impl(logx, Sum);
    }else if(type=="col"){
    column_sum_impl(logx, Sum); 
    }
  }
    
  if(log==0){
    if(type=="row"){
    row_sum_impl(x, Sum);
    }else if(type=="col"){
    column_sum_impl(x, Sum);  
    }
  }
  
  
}


template<typename T> 
void rowsum_Templated(
    Rcpp::S4  xR, 
    Rcpp::S4  SumR,
    std::string type,
    int log){    //   Rcpp::IntegerVector numWorkItems){
  
  const bool BisVCL=1;
  const int ctx_id = INTEGER(xR.slot(".context_index"))[0]-1;
  
  
  std::shared_ptr<viennacl::matrix<T> > x =getVCLptr<T>(xR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::vector_base<T> > Sum = getVCLVecptr<T>(SumR.slot("address"), BisVCL, ctx_id);
  rowsum<T>(*x,*Sum,type,log);
  
  //return (Rcpp::wrap(1L));
}




// [[Rcpp::export]]
void rowsumBackend(
    Rcpp::S4  xR, 
    Rcpp::S4  SumR,
    std::string type,
    int log){      //  Rcpp::IntegerVector numWorkItems){
  
  Rcpp::traits::input_parameter< std::string >::type classVarR(RCPP_GET_CLASS(xR));
  std::string precision_type = (std::string) classVarR;
  
  if(precision_type == "fvclMatrix") {
    return (rowsum_Templated<float>(xR,SumR, type,log));
  } else if (precision_type == "dvclMatrix") {
    return (rowsum_Templated<double>(xR, SumR, type, log));
  } else {
    Rcpp::warning("class of var must be fvclMatrix or dvclMatrix");
  }

}






