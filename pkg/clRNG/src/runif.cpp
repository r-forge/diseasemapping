#include "clRNG.hpp"
#include <clRNG/mrg31k3p.h>
#include <CL/mrg31k3pkernelStringSeparate.hpp>
#include <CL/mrg31k3pkernelStringCommon.hpp>

using namespace Rcpp;
using namespace viennacl;	
using namespace viennacl::linalg;

// clRNG -> Matrix
void convertclRngMat(clrngMrg31k3pStream* streams, Rcpp::IntegerMatrix result) {
  
  int numWorkItems = result.nrow();
  int Ditem,Delement,Dcis,Dg;
  for(Ditem =0;Ditem < numWorkItems;Ditem++){
    for(Delement=0;Delement < 3;Delement++){
      
      Dcis=0;
      Dg=0;
      result(Ditem, Dcis*6 + Dg*3 + Delement) = streams[Ditem].current.g1[Delement];//0,0; 0,1
      Dg=1;
      result(Ditem,Dcis*6 + Dg*3 + Delement) = streams[Ditem].current.g2[Delement];//0,3; 0,4
      
      Dcis=1;
      Dg=0;
      result(Ditem,Dcis*6 + Dg*3 + Delement) = streams[Ditem].initial.g1[Delement];//0,6; 0,7
      Dg=1;
      result(Ditem,Dcis*6 + Dg*3 + Delement) = streams[Ditem].initial.g2[Delement];//0,9; 0,10
      
      Dcis=2;
      Dg=0;
      result(Ditem,Dcis*6 + Dg*3 + Delement) = streams[Ditem].substream.g1[Delement];//0,12
      Dg=1;
      result(Ditem,Dcis*6 + Dg*3 + Delement) = streams[Ditem].substream.g2[Delement];//0,15
      
    }
  }
  
}

//matrix ->clRNG streams
void convertMatclRng(Rcpp::IntegerMatrix Sin, clrngMrg31k3pStream* streams){
  
  int Ditem,Delement,Dcis,Dg;
  int numWorkItems = Sin.ncol();
  
  for(Ditem =0;Ditem < numWorkItems;Ditem++){
    for(Delement=0;Delement < 3;Delement++){
      
      Dcis=0;
      Dg=0;
      streams[Ditem].current.g1[Delement] = Sin(Ditem,Dcis*6 + Dg*3 + Delement);
      Dg=1;
      streams[Ditem].current.g2[Delement] = Sin(Ditem,Dcis*6 + Dg*3 + Delement);
      
      Dcis=1;
      Dg=0;
      streams[Ditem].initial.g1[Delement] = Sin(Ditem,Dcis*6 + Dg*3 + Delement);
      Dg=1;
      streams[Ditem].initial.g2[Delement]=Sin(Ditem,Dcis*6 + Dg*3 + Delement);
      
      Dcis=2;
      Dg=0;
      streams[Ditem].substream.g1[Delement]=Sin(Ditem,Dcis*6 + Dg*3 + Delement);
      Dg=1;
      streams[Ditem].substream.g2[Delement] = Sin(Ditem,Dcis*6 + Dg*3 + Delement);
    }
  }
  
  
}


template <typename T> 
std::string mrg31k3pTypeString() {
  return("undefined");}

template <> std::string mrg31k3pTypeString<double>(){
  std::string result;
  result = 
    "\n#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n\n" + 
    mrg31k3pCommon +
    mrg31k3pTemplateStringFirst + "double" +
    mrg31k3pTemplateStringSecond + 
  "4.656612873077392578125e-10 * " +  /* 1/2^31 */
  mrg31k3pTemplateStringThird;
  return(result);
}

template <> std::string mrg31k3pTypeString<float>(){
  std::string result;
  result =  mrg31k3pCommon + 
    mrg31k3pTemplateStringFirst + "float" +
    mrg31k3pTemplateStringSecond + 
    "4.6566126e-10 * " +  /* 1/2^31 */
    mrg31k3pTemplateStringThird;
  return(result);
}


template <> std::string mrg31k3pTypeString<int>(){
  std::string result;
  result = mrg31k3pCommon + 
    mrg31k3pTemplateStringFirst + "int" +
    mrg31k3pTemplateStringSecond + "(int) " + 
    mrg31k3pTemplateStringThird;
  ;
  return(result);
}

//////////////////////////////////the main function;//////////////////////
template<typename T>
void runifGpu(
    viennacl::vector_base<T> &x, 
    Rcpp::IntegerMatrix streamsR, 
    int numWorkItems,
    int numLocalItems,
    int ctx_id){
  
  
  viennacl::ocl::context ctx(viennacl::ocl::get_context(ctx_id));
  std::string mrg31k3pkernelString = mrg31k3pTypeString<T>();

//  Rcout << "hi\n\n" <<mrg31k3pkernelString << "\n\n";

  // create streams
  size_t streamBufferSize;   
  clrngStatus err;
  
  //Reserve memory space for count stream objects, without creating the stream objects. Returns a pointer to the newly allocated buffer. 
  clrngMrg31k3pStream* streams = clrngMrg31k3pAllocStreams(numWorkItems, &streamBufferSize, &err);
  //count	Number of stream objects to allocate.
  
  
  
  // transfer streams to opencl as matrix
  // convert to crngMgr31k3pStream in opencl, but still on host
  convertMatclRng(streamsR, streams);
  
  
  // Create buffer to transfer streams to the device.
  viennacl::vector<char> bufIn(ctx.create_memory_without_smart_handle( 
      CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,  streamBufferSize, (void *) streams), 1);
  
  
  // kernel, in the kernel will Copy RNG host stream objects from global memory into private memory
  viennacl::ocl::program &my_prog = ctx.add_program(mrg31k3pkernelString, "my_kernel");

  viennacl::ocl::kernel &randomUniform = my_prog.get_kernel("mrg31k3p");
  randomUniform.global_work_size(0, numWorkItems);
  randomUniform.local_work_size(0, numLocalItems);
  
  int Nsim = x.size();
  viennacl::ocl::enqueue(randomUniform(bufIn, x, Nsim) ); //streams, out, vector_size
  
  
  // copy streams back to cpu
  viennacl::backend::memory_read(bufIn.handle(), 0, streamBufferSize,streams);

    // then transfer to R object, //return streams to R 
  convertclRngMat(streams, streamsR);
  
}

////////////////////////////////////////////////////////////////////////////////////////////////
void runifGpuHost(viennacl::vector_base<double> &x)//use them to generate numbers on the host
{
  
  int D, N=x.size();
  clrngMrg31k3pStream* stream = clrngMrg31k3pCreateStreams(NULL, 1, NULL, NULL);
  
  for(D =0; D<N; D++) {
    x(D) = clrngMrg31k3pRandomU01(stream);
  }
}
////////////////////////////////////////////////////////////////////////////////////////////////

//[[Rcpp::export]]
Rcpp::IntegerMatrix  cpp_mrg31k3pCreateStreams(int numWorkItems) 
{
  
  Rcpp::IntegerMatrix result=Rcpp::IntegerMatrix(numWorkItems,18L);
  
  colnames(result) = CharacterVector::create(
    "current.g1.1", "current.g1.2", "current.g1.3", "current.g2.1", "current.g2.2", "current.g2.3",
    "initial.g1.1", "initial.g1.2", "initial.g1.3", "initial.g2.1", "initial.g2.2", "initial.g2.3",
    "substream.g1.1", "substream.g1.2", "substream.g1.3", "substream.g2.1", "substream.g2.2", "substream.g2.3");
  
  size_t streamBufferSize;
  clrngStatus err;
  
  int Ditem,Delement,Dcis,Dg;
  
  clrngMrg31k3pStream* streams = clrngMrg31k3pCreateStreams(NULL, numWorkItems, &streamBufferSize, &err);//line 299 in mrg31k3p.c
  
  convertclRngMat(streams, result);
  
  return result;
}





/////////////////////////////////uniforms on GPU function;//////////
template<typename T> 
SEXP runifGpu(
    Rcpp::S4  xR,   //vector
    Rcpp::IntegerMatrix streamsR,   //vector
    int max_global_size,     
    int max_local_size) 
{
  // data
  const bool BisVCL=1;
  const int ctx_id = INTEGER(xR.slot(".context_index"))[0]-1;
  
  std::shared_ptr<viennacl::vector_base<T> > x = getVCLVecptr<T>(xR.slot("address"), BisVCL, ctx_id);
  
  runifGpu<T>(*x, streamsR, max_global_size, max_local_size, ctx_id);
  
  return(Rcpp::wrap(1L));	
}

//[[Rcpp::export]]
SEXP cpp_runifGpu(
    Rcpp::S4  xR,   //vector
    Rcpp::IntegerMatrix streamsR,   //vector
    int max_global_size,     
    int max_local_size,
    std::string type) 
{
  if(type == "float") {
    return runifGpu<float>(xR, streamsR, max_global_size, max_local_size);
  } else if (type=="integer") {
    return runifGpu<int>(xR, streamsR, max_global_size, max_local_size);
  } else {
    return runifGpu<double>(xR, streamsR, max_global_size, max_local_size);
  }
}



