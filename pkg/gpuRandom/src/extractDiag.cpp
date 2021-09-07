#include "gpuRandom.hpp"
#define DEBUGKERNEL
template <typename T> 
std::string extract_some_diag_string(int Ndatasets,    // number of Y's
                                     int NpadColY, 
                                     int NpadColYX,
                                     int NpadBetweenMatricesYX) {  //ssqYX.size2()
  
  std::string typeString = openclTypeString<T>();
  
  
  std::string result = "";
  
  if(typeString == "double") {
    result += "\n#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n";
  }
  
  
  result += 
    //"#define NthisIteration " + std::to_string(NthisIteration) + "\n"
    "#define Ndatasets " + std::to_string(Ndatasets) + "\n"    
    "#define NpadColY " + std::to_string(NpadColY) + "\n"
    "#define NpadBetweenMatricesYX "   + std::to_string(NpadBetweenMatricesYX) + "\n"
    "#define NpadColYX "   + std::to_string(NpadColYX) + "\n";
  
  result += 
    "\n__kernel void extract_some_diag(\n"
    "__global " + typeString + " *ssqY,\n"
    "__global " + typeString + " *ssqYX,\n"  
    " int Nrowstart,\n"
    " int NthisIteration ){\n\n";
  
  result +=
    "int Drow, Dcol;\n"
    "int NrowEnd = Nrowstart+NthisIteration;\n";
  
  result += 
    
    "for (Drow = Nrowstart+get_global_id(0); Drow < NrowEnd; Drow += get_global_size(0)) {\n"
    "for (Dcol = get_global_id(1); Dcol < Ndatasets; Dcol += get_global_size(1)) {\n" 
    
    "ssqY[Drow * NpadColY+Dcol] = ssqYX[ Drow*NpadBetweenMatricesYX + Dcol* NpadColYX + Dcol ];\n"
    
    "}\n"
    "}\n"
    "}\n";
  
  return(result);
}











void extractDiag(
    viennacl::matrix<double> &ssqYX,// must be batches of square matrices
    viennacl::matrix<double> &ssqY,
    Rcpp::IntegerVector NparamPerIter,
    int Ndatasets,
    Rcpp::IntegerVector numWorkItems,
    int ctx_id) {
  

  
  std::string extractSomeDiagKernelString = extract_some_diag_string<double>(//NparamPerIter[0],  
                 Ndatasets,
                 ssqY.internal_size2(), 
                 ssqYX.internal_size2(), 
                 ssqYX.internal_size2()*ssqYX.size2());

  
#ifdef DEBUGKERNEL
  Rcpp::Rcout << extractSomeDiagKernelString << "\n\n";
#endif    
  
  
  
  // the context
  viennacl::ocl::context ctx(viennacl::ocl::get_context(ctx_id));
  viennacl::ocl::program & my_prog = ctx.add_program(extractSomeDiagKernelString, "my_kernel");

  
  viennacl::ocl::kernel &extractSomeDiagKernel = my_prog.get_kernel("extract_some_diag");
  extractSomeDiagKernel.global_work_size(0, numWorkItems[0]);
  extractSomeDiagKernel.global_work_size(1, numWorkItems[1]);
  
  extractSomeDiagKernel.local_work_size(0, 1L);
  extractSomeDiagKernel.local_work_size(1, 1L);
  
  int Nparams = ssqYX.size1()/ssqYX.size2();
  int Niter = ceil( ( (double) Nparams) / ((double) NparamPerIter[0]));
  int Diter;
  int endThisIteration;
  int DiterIndex, NthisIteration;
  
  
  for (Diter=0, DiterIndex=0; 
       Diter< Niter; 
       Diter++, DiterIndex += NparamPerIter[0]){
    
    endThisIteration = std::min(DiterIndex + NparamPerIter[0], Nparams);
    NthisIteration = endThisIteration - DiterIndex;
    
    // if(verbose[0]) {
    //   Rcpp::Rcout << "\n" << "Diter " << Diter <<" DiterIndex " << DiterIndex << " endThisIteration " << 
    //     endThisIteration << " Nthisiteration " << NthisIteration  << "\n";
    // }
    
    
    viennacl::ocl::enqueue(extractSomeDiagKernel(ssqY, ssqYX, DiterIndex, NthisIteration));     //Rcpp::Rcout << "after enqueue" << "\n\n";

    
  } // Diter
  
  
  
  
  
  
}



//template<typename T> 
void extractDiagTemplated(
    Rcpp::S4 ssqYXR,// viennacl::vector_base<int>  rowSum, viennacl::vector_base<int>  colSum,  
    Rcpp::S4 ssqYR,
    Rcpp::IntegerVector NparamPerIter,
    int Ndatasets,
    Rcpp::IntegerVector numWorkItems) {

  
  const bool BisVCL=1;
  const int ctx_id = INTEGER(ssqYXR.slot(".context_index"))[0]-1;
  std::shared_ptr<viennacl::matrix<double> > ssqYX = getVCLptr<double>(ssqYXR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<double> > ssqY = getVCLptr<double>(ssqYR.slot("address"), BisVCL, ctx_id);
  
  extractDiag(
    *ssqYX,// viennacl::vector_base<int>  rowSum, viennacl::vector_base<int>  colSum,  
    *ssqY,
    NparamPerIter,
    Ndatasets,
    numWorkItems,
    ctx_id);

  
}






// [[Rcpp::export]]
void extractDiagBackend(
    Rcpp::S4 ssqYXR,// viennacl::vector_base<int>  rowSum, viennacl::vector_base<int>  colSum,  
    Rcpp::S4 ssqYR,
    Rcpp::IntegerVector NparamPerIter,
    int Ndatasets,
    Rcpp::IntegerVector numWorkItems) {
  
  
  extractDiagTemplated(ssqYXR,ssqYR, NparamPerIter,Ndatasets,numWorkItems);
  
  
}
