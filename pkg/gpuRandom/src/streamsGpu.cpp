#include "gpuRandom.hpp"
#define DEBUGKERNEL

std::string streamsString(int NpadStreams) {  
  
  
  std::string result = "";
  
  
  result += 
    "#define NpadStreams " + std::to_string(NpadStreams) + "\n"
    "#define mrg31k3p_M1 2147483647 \n"
    "#define mrg31k3p_M2 2147462579 \n";
  
  result += 
    "\n__kernel void createStreams(\n" 
    "__global int *creatorInitialGlobal, \n"
    "__global int *streams,\n"
    "int Nstreams){\n";
  
  
  result +=
    "int creatorNextState[6], creatorCurrentState[6], g[6];\n"  
    "int i, row, col, Dstream;\n"
    "float acc; \n"
    "const int jmatrix[18]={1702500920, 1849582496, 1656874625,\n"
    " 828554832, 1702500920, 1512419905,\n"
    " 1143731069,  828554832,  102237247,\n"
    " 796789021, 1464208080,  607337906, \n"
    " 1241679051, 1431130166, 1464208080, \n"
    " 1401213391, 1178684362, 1431130166};\n\n\n";
  
  
  result += 
    
    // copy creatorInitialGlobal to creatorCurrentState
    " for (i=0; i<6; i++) {\n"
    "creatorNextState[i] = creatorInitialGlobal[i];\n"
    "}\n"
    
    "for(Dstream = 0;   Dstream < Nstreams;    Dstream++){\n"
    
    // upate creatorNext from creatorCurrentState,
    "for (i=0; i<6; i++) {\n"
    " streams[Dstream * NpadStreams +  i] = \n"//initial
    "streams[Dstream * NpadStreams + 6 + i] = \n"//current
    "streams[Dstream * NpadStreams + 12 + i] = \n"// substream
    "creatorCurrentState[i] = creatorNextState[i];\n"
    "}\n"
    
    
    
    // Matrix-vector modular multiplication
    // modMatVec(creator->nuA1, creator->nextState.g1, creator->nextState.g1, mrg31k3p_M1);
    "for (row=0; row<3; row++){\n"
    "acc = 0.0;\n"
    "for (col=0; col<3; col++){\n"
    "acc += jmatrix[3 * row + col] * creatorNextState[col];\n"
    " }\n"
    "creatorNextState[row] = fmod((float)acc, (float)mrg31k3p_M1);\n"
    "}\n"
    
    // modMatVec(creator->nuA2, creator->nextState.g2, creator->nextState.g2, mrg31k3p_M2);
    "for (row=3; row<6; row++){\n"
    " acc = 0.0;\n"
    "for (col=0; col<3; col++){\n"
    "acc += jmatrix[3 * row + col] * creatorNextState[col+3];\n"
    "}\n"
    "creatorNextState[row] = fmod((float)acc, (float)mrg31k3p_M2);\n"
    "}\n"
    
    "}\n" // loop through streams
    
    
    // copy current state to substream
    
    
    
    "}\n"; 
  
  return(result);
}








void CreateStreamsGpu(
    viennacl::matrix_base<int> &streams,
    int ctx_id) {
  
  std::vector<int>   cpu_vector = {12345, 12345, 12345, 12345, 12345, 12345};
  viennacl::vector_base<int> creatorInitialGlobal(6);
  /* fill a vector on CPU
   for (size_t i=0; i<6; ++i)
   cpu_vector[i] = 12345; */
  
  // fill a vector on GPU with data from CPU - faster versions:
  copy(cpu_vector, creatorInitialGlobal);  //option 1 // copy(cpu_vector.begin(), cpu_vector.end(), vcl_vector.begin()); //option 2
  
  const int Nstreams = streams.size1(); //# of rows
  std::string streamsKernelString = streamsString(
    streams.internal_size2() 
  );
  
  // the context
  viennacl::ocl::switch_context(ctx_id);
  viennacl::ocl::program & my_prog = viennacl::ocl::current_context().add_program(streamsKernelString, "my_kernel");
  
#ifdef DEBUGKERNEL
  Rcpp::Rcout << streamsKernelString << "\n\n";
#endif  
  
  viennacl::ocl::kernel &streamsKernel = my_prog.get_kernel("createStreams");
  streamsKernel.global_work_size(0, 1L);
  streamsKernel.global_work_size(1, 1L);
  
  //streamsKernel.local_work_size(0, 1L);
  //streamsKernel.local_work_size(1, 1L);
  Rcpp::Rcout << "11" << "\n\n";
  
  viennacl::ocl::enqueue(streamsKernel(creatorInitialGlobal, streams, Nstreams) );
  Rcpp::Rcout << "22" << "\n\n";
}



SEXP CreateStreamsGpuTemplated(
    Rcpp::S4 streamsR){
  
  const bool BisVCL=1;
  const int ctx_id = INTEGER(streamsR.slot(".context_index"))[0]-1;
  std::shared_ptr<viennacl::matrix_base<int> > streams = getVCLptr<int>(streamsR.slot("address"), BisVCL, ctx_id);
  std::cout<< "33\n";    
  CreateStreamsGpu(*streams, ctx_id);
  
  return Rcpp::wrap(0L);
  std::cout<< "44\n";  
}





//[[Rcpp::export]]
SEXP CreateStreamsGpuBackend(
    Rcpp::S4 streamsR) {
  
  
  Rcpp::traits::input_parameter< std::string >::type classVarR(RCPP_GET_CLASS(streamsR));
  std::string precision_type = (std::string) classVarR;
  
  if(precision_type == "fvclMatrix" || precision_type == "dvclMatrix") {
    Rcpp::warning("must be matrix of integers!\n");
  }
  Rcpp::Rcout << "55" << "\n\n";
  if(precision_type == "ivclMatrix") {
    CreateStreamsGpuTemplated(streamsR);
  } 
  
  
}





