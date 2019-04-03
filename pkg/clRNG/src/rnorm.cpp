#include "clRNG.hpp"
#include <clRNG/mrg31k3p.h>
#include <CL/mrg31k3pkernelString.hpp>

using namespace Rcpp;
using namespace viennacl; 
using namespace viennacl::linalg;



template <typename T> 
std::string mrg31k3pStringType() {
	return("undefined");}
template <> std::string mrg31k3pStringType<double>(){
	return("mrg31k3pDoubleNorm");}
template <> std::string mrg31k3pStringType<float>(){
	return("mrg31k3pFloatNorm");}



template<typename T>
void rnormGpu(
	viennacl::vector_base<T> &x, 
  Rcpp::IntegerMatrix streamsR, 
	int numWorkItems,    
	int numLocalItems,
	int ctx_id){

viennacl::ocl::context ctx(viennacl::ocl::get_context(ctx_id));

// create streams
size_t streamBufferSize;   
clrngStatus err;

//Reserve memory space for count stream objects, without creating the stream objects. Returns a pointer to the newly allocated buffer. 
clrngMrg31k3pStream* streams = clrngMrg31k3pAllocStreams(numWorkItems, &streamBufferSize, &err);//count	Number of stream objects to allocate.


// transfer streams to opencl as matrix
// convert to crngMgr31k3pStream in opencl, but still on host
convertMatclRng(streamsR, streams);

// Create buffer to transfer streams to the device.
viennacl::vector<char> bufIn(ctx.create_memory_without_smart_handle(
CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,  streamBufferSize, (void *) streams), 1);


// kernel, in the kernel will Copy RNG host stream objects from global memory into private memory
viennacl::ocl::program &my_prog = ctx.add_program(mrg31k3pkernelString, "my_kernel");
viennacl::ocl::kernel &randomNorm = my_prog.get_kernel(mrg31k3pStringType<T>());

randomNorm.global_work_size(0, numWorkItems[0]);
randomNorm.global_work_size(1, numWorkItems[1]);

randomNorm.local_work_size(0, numLocalItems[0]);
randomNorm.local_work_size(1, numLocalItems[1]);


int Nsim = x.size();
viennacl::ocl::enqueue(randomNorm(bufIn, x, Nsim) ); //streams, out, vector_size

// copy streams back to cpu
viennacl::backend::memory_read(bufIn.handle(), 0, streamBufferSize,streams);

// then transfer to R object, //return streams to R 
convertclRngMat(streams, streamsR);
}

/////////////////////////////////random normals on GPU function;//////////
template<typename T> 
SEXP rnormGpu(
    Rcpp::S4  xR,   //vector
    Rcpp::IntegerMatrix streamsR,   //vector
    int max_global_size,     
    int max_local_size) 
{
  // data
  const bool BisVCL=1;
  const int ctx_id = INTEGER(xR.slot(".context_index"))[0]-1;
  
  std::shared_ptr<viennacl::vector_base<T> > x = getVCLVecptr<T>(xR.slot("address"), BisVCL, ctx_id);
  
  rnormGpu<T>(*x, streamsR, max_global_size, max_local_size, ctx_id);
  
  return(Rcpp::wrap(1L));	
}

//[[Rcpp::export]]
SEXP cpp_rnormGpu(
    Rcpp::S4  xR,   //vector
    Rcpp::IntegerMatrix streamsR,   //vector
    int max_global_size,     
    int max_local_size,
    std::string type) 
{
  if(type == "float") {
    return rnormGpu<float>(xR, streamsR, max_global_size, max_local_size);
  } else if (type=="double") {
    return rnormGpu<double>(xR, streamsR, max_global_size, max_local_size);
  } 
}



















