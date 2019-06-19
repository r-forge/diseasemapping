#include <CL/fisher_sim.hpp>   
#include "random_number.hpp"

using namespace Rcpp;
using namespace viennacl;	
using namespace viennacl::linalg;



void fisher_sim_gpu(
    viennacl::vector_base<int> &sr, 
    viennacl::vector_base<int> &sc,
    viennacl::vector_base<double> &ans, 
    Rcpp::IntegerMatrix streamsR, 
    Rcpp::IntegerVector numWorkItems,
    int ctx_id){
  
  
  viennacl::ocl::context ctx(viennacl::ocl::get_context(ctx_id));
  std::string kernel_string= FisherSimkernelString;
  
  //Rcout << FisherSimkernelString << "\n\n";
  // create streams
  size_t streamBufferSize;   
  clrngStatus err;
  
  
  //Reserve memory space for count stream objects, without creating the stream objects. 
  //Returns a pointer to the newly allocated buffer. 
  clrngMrg31k3pStream* streams = clrngMrg31k3pAllocStreams(numWorkItems[0]*numWorkItems[1], &streamBufferSize, &err);
  //count	Number of stream objects to allocate.
  //#ifdef UNDEF   
  // transfer streams to opencl as matrix
  // convert to crngMgr31k3pStream in opencl, but still on host
  convertMatclRng(streamsR, streams);
  
  // Create buffer to transfer streams to the device.
  viennacl::vector<char> bufIn(ctx.create_memory_without_smart_handle( 
      CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,  streamBufferSize, (void *) streams), 1);
  
  // kernel, in the kernel will Copy RNG host stream objects from global memory into private memory
  // add kernel to program

 /* std::string stuff = "\n\n__kernel void fisher_sim_gpu2(\n"
  "int *nrow\n"
  ") { \n"
  "   double ans;\n"
  "};\n";*/


  viennacl::ocl::program &my_prog = ctx.add_program(kernel_string, "my_kernel");
  // get compiled kernel function
  viennacl::ocl::kernel &fisher_sim = my_prog.get_kernel("fisher_sim_gpu");
  
  fisher_sim.global_work_size(0, numWorkItems[0]);
  fisher_sim.global_work_size(1, numWorkItems[1]);
  
  
  const int nr = sr.size(), nc = sc.size();

  scalar<int> nScalar;
  //viennacl::linalg::sum_impl(sr, &n);  //User interface function for computing the sum of all elements of a vector. 
  sum_impl(sr, nScalar);
  const int n = nScalar;
  int vsize= ans.size();
  
  viennacl::vector_base<int> observed = viennacl::vector_base<int>(nr*nc, ctx); 

  viennacl::vector_base<double> fact = viennacl::vector_base<double>(n+1, ctx); 
  // TO DO: put this in local memory
  viennacl::vector_base<int> jwork = viennacl::vector_base<int>(nc, ctx); 


  viennacl::ocl::enqueue(fisher_sim(
    nr, nc, 
    sr, sc, 
    n, vsize, 
    observed, fact, jwork, 
    ans,
    bufIn
    ) ); //streams, out, vector_size



  
  // copy streams back to cpu
  viennacl::backend::memory_read(bufIn.handle(), 0, streamBufferSize, streams);
  //#endif 
  // then transfer to R object, //return streams to R 
  convertclRngMat(streams, streamsR);
}



//[[Rcpp::export]]
SEXP cpp_fisher_sim_gpu(
    Rcpp::S4  srR, 
    Rcpp::S4  scR,
    Rcpp::S4  ansR,
    Rcpp::IntegerMatrix streamsR,   //vector
    Rcpp::IntegerVector max_global_size){
  
  //const bool BisVCL=1;
  const int ctx_id = INTEGER(srR.slot(".context_index"))[0]-1;
  std::shared_ptr<viennacl::vector_base<int> > sr = getVCLVecptr<int>(srR.slot("address"), 1, ctx_id);
  std::shared_ptr<viennacl::vector_base<int> > sc = getVCLVecptr<int>(scR.slot("address"), 1, ctx_id);
  std::shared_ptr<viennacl::vector_base<double> > ans = getVCLVecptr<double>(ansR.slot("address"), 1, ctx_id);
  
  
  fisher_sim_gpu(*sr, *sc, *ans, streamsR, max_global_size, ctx_id);
  
  return (ansR);
}
