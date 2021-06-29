#include "gpuRandom.hpp"
#define DEBUGKERNEL

std::string streamsString(int NpadStreams, 
                          const int keepinitial) {  
  
  
  std::string result = "";
  
  
  result += 
    "#define NpadStreams " + std::to_string(NpadStreams) + "\n"
    "#define mrg31k3p_M1 2147483647 \n"
    "#define mrg31k3p_M2 2147462579 \n\n\n";
  

   result += "__constant uint jmatrix[18]= {1702500920, 1849582496, 1656874625,\n"
   " 828554832, 1702500920, 1512419905,\n"
   " 1143731069,  828554832,  102237247,\n"
   " 796789021, 1464208080,  607337906, \n"
   " 1241679051, 1431130166, 1464208080, \n"
   " 1401213391, 1178684362, 1431130166};\n\n\n";
 
 
 /* 
  result += "__constant cl_uint jmatrix[18]= {1, 2, 3,\n"
  " 4, 5, 6,\n"
  " 7,  8, 9,\n"
  "1, 0, 0,\n"
  " 0, 1, 0,\n"
  " 0,  0,  1 };\n\n\n";
  */  
  
  
  result += 
    "\n__kernel void createStreams(\n"    
    "__global uint *creatorInitialGlobal, \n"
    "__global uint *streams,\n"
    "int Nstreams){\n";
  
  
  result +=
    "uint creatorNextState[6], creatorCurrentState[6];\n"  
    "int i, row, col, Dstream;\n"
    "uint acc; \n";
  
  result +=  
    " for (i=0; i<6; i++) {\n"
    "creatorCurrentState[i] = creatorInitialGlobal[i];\n"
    "}\n";
  
  
  
  
  /*
   if (keepinitial == 0) {
   
   result +=     
   // Matrix-vector modular multiplication
   // modMatVec(creator->nuA1, creator->nextState.g1, creator->nextState.g1, mrg31k3p_M1);
   "for (row=0; row<3; row++){\n"
   "acc = 0;\n"
   "for (col=0; col<3; col++){\n"
   "acc += jmatrix[3 * row + col] * creatorInitialGlobal[col];\n"
   " }\n"
   //"creatorCurrentState[row] = acc; \n"
   "creatorCurrentState[row] = acc % mrg31k3p_M1;\n"
   // "creatorCurrentState[row] = fmod((double)acc, (double)mrg31k3p_M1);\n"
   "}\n"
   
   
   // modMatVec(creator->nuA2, creator->nextState.g2, creator->nextState.g2, mrg31k3p_M2);
   "for (row=3; row<6; row++){\n"
   " acc = 0;\n"
   "for (col=0; col<3; col++){\n"
   "acc += jmatrix[3 * row + col] * creatorInitialGlobal[col+3];\n"
   "}\n"
   //"creatorCurrentState[row] = acc; \n"
   "creatorCurrentState[row] = acc % mrg31k3p_M2;\n"
   // "creatorCurrentState[row] = fmod((float)acc, (float)mrg31k3p_M2);\n"
   "}\n\n";
   
   }
   
   */ 
  
  
  
  
  result +=  
  
  "for(Dstream = 0;   Dstream < Nstreams;    Dstream++){\n\n"
  
  // upate creatorNext from creatorCurrentState,
  "for (i=0; i<6; i++) {\n"
  " streams[Dstream * NpadStreams +  i] = \n"//initial
  "streams[Dstream * NpadStreams + 6 + i] = \n"//current
  "streams[Dstream * NpadStreams + 12 + i] = "// substream
  " creatorNextState[i] = creatorCurrentState[i];\n"
  "}\n"
  
  
  
  // Matrix-vector modular multiplication
  // modMatVec(creator->nuA1, creator->nextState.g1, creator->nextState.g1, mrg31k3p_M1);
  "for (row=0; row<3; row++){\n"
  "acc = 0;\n"
  "for (col=0; col<3; col++){\n"
  "acc += jmatrix[3 * row + col] * creatorNextState[col];\n"
  " }\n"
  "creatorCurrentState[row] = acc; \n"
  // "creatorCurrentState[row] = acc % mrg31k3p_M1;\n"
  // "creatorCurrentState[row] = fmod((float)acc, (float)mrg31k3p_M1);\n"
  "}\n"
  
  
  
  // modMatVec(creator->nuA2, creator->nextState.g2, creator->nextState.g2, mrg31k3p_M2);
  "for (row=3; row<6; row++){\n"
  " acc = 0;\n"
  "for (col=0; col<3; col++){\n"
  "acc += jmatrix[3 * row + col] * creatorNextState[col+3];\n"
  "}\n"
   "creatorCurrentState[row] = acc; \n"
   //"creatorCurrentState[row] = acc % mrg31k3p_M2;\n"
  // "creatorCurrentState[row] = fmod((float)acc, (float)mrg31k3p_M2);\n"
  "}\n"
  
  "}\n" // loop through streams
  
  "}\n"; 
  
  return(result);
}








void CreateStreamsGpu(
    viennacl::vector_base<cl_uint> &creatorInitialGlobal,
    viennacl::matrix_base<cl_uint> &streams, 
    const int keepinitial,
    int ctx_id) {
  
  /* std::vector<int>   cpu_vector = {12345, 12345, 12345, 12345, 12345, 12345};
   viennacl::vector_base<int> creatorInitialGlobal(6); */
  /* fill a vector on CPU
   for (size_t i=0; i<6; ++i)
   cpu_vector[i] = 12345; */
  
  // fill a vector on GPU with data from CPU - faster versions:
  //copy(cpu_vector, creatorInitialGlobal);  //option 1 // copy(cpu_vector.begin(), cpu_vector.end(), vcl_vector.begin()); //option 2
  
  const int Nstreams = streams.size1(); //# of rows
  std::string streamsKernelString = streamsString(
    streams.internal_size2(), 
    keepinitial
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
  //Rcpp::Rcout << "11" << "\n\n";
  
  viennacl::ocl::enqueue(streamsKernel(creatorInitialGlobal, streams, Nstreams) );
  
  Rcpp::Rcout << streams(0,0) << "\n" << streams(0,1) << "\n"<< streams(0,2) << "\n\n";
  Rcpp::Rcout << streams(1,0) << "\n" << streams(1,1) << "\n"<< streams(1,2) << "\n\n";
  Rcpp::Rcout << streams(2,0) << "\n" << streams(2,1) << "\n"<< streams(2,2) << "\n\n";
}










void CreateStreamsGpuTemplated(
    Rcpp::S4 creatorInitialGlobalR,
    Rcpp::S4 streamsR,
    const int keepinitial){
  
  const bool BisVCL=1;
  const int ctx_id = INTEGER(streamsR.slot(".context_index"))[0]-1;
  std::shared_ptr<viennacl::vector_base<cl_uint> > creatorInitialGlobal = getVCLVecptr<cl_uint>(creatorInitialGlobalR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix_base<cl_uint> > streams = getVCLptr<cl_uint>(streamsR.slot("address"), BisVCL, ctx_id);
  
  
  
  //std::cout<< "33\n";    
  CreateStreamsGpu(*creatorInitialGlobal, *streams, keepinitial, ctx_id);
  
  //std::cout<< "44\n";  
  //return Rcpp::wrap(0L);
}





//[[Rcpp::export]]
void CreateStreamsGpuBackend(
    Rcpp::S4 creatorInitialGlobalR,    
    Rcpp::S4 streamsR,
    const int keepinitial){
  
  
  /*
   Rcpp::traits::input_parameter< std::string >::type classVarR(RCPP_GET_CLASS(creatorInitialGlobalR));
   std::string precision_type = (std::string) classVarR;
   */
  
  //Rcpp::traits::input_parameter< std::string >::type classstream(RCPP_GET_CLASS(streamsR));
  //std::string precision_type_stream = (std::string) classstream;
  
  
  //Rcpp::Rcout << "55" << "\n\n";
  
  
  
  CreateStreamsGpuTemplated(creatorInitialGlobalR, streamsR, keepinitial);
  
}










