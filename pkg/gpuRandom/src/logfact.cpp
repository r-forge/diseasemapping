#include "gpuRandom.hpp"


#define DEBUGKERNEL

template <typename T> 
std::string logfactString() {  
  

  std::string typeString = openclTypeString<T>();  // type of the log factorial
  
  std::string result = "";
  
  if(typeString == "double") {
    
    result += "\n#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n";
    
  }
  

  
   result += 
   "\n\n__kernel void logfactorial(\n"
//   "  __global int *vector,\n"  
   "  __global " + typeString + "* out,\n"  
   " const int numElements\n"
   "){\n";  
   
   
   
   result += "const int size = get_global_size(1)*get_global_size(0);\n"
   "const int index = get_global_id(1)*get_global_size(0) + get_global_id(0);\n"
   "int D;\n";
    result += typeString + " Dfloat;\n\n";   
   
   result += "for(D=index; D< numElements; D+=size){\n"
    " Dfloat = D;\n"
     "out[D]=lgamma(Dfloat);\n"
   "}\n"
   "}\n";

  return(result);
}


//template <typename T> 
void logfactorial(){
//    viennacl::vector<T>  &output, //must be same length with x
    Rcpp::IntegerVector numWorkItems,
    int ctx_id) {
#ifdef UNDEF  
  const int numelements=output.size();
  
  std::string logKernelString = logfactString<T>();
  
  // the context
  viennacl::ocl::switch_context(ctx_id);
  viennacl::ocl::program & my_prog = viennacl::ocl::current_context().add_program(logKernelString, "my_kernel");
  
#ifdef DEBUGKERNEL
  Rcpp::Rcout << logKernelString << "\n\n";
#endif  
  
  
  viennacl::ocl::kernel &lfactorialKernel = my_prog.get_kernel("logfactorial");
  
  lfactorialKernel.global_work_size(0, numWorkItems[0]);
  lfactorialKernel.global_work_size(1, numWorkItems[1]);
  
 // lfactorialKernel.local_work_size(0, 1L);
 // lfactorialKernel.local_work_size(1, 1L);
  
  
  viennacl::ocl::enqueue(lfactorialKernel(output, numelements) );
#endif  
  
}

















