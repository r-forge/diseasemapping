#include <CL/fisher_sim.hpp>   
#include "random_number.hpp"

#define DEBUG

using namespace Rcpp;
using namespace viennacl;	
using namespace viennacl::linalg;


#include <string>

std::string FisherSimkernelString(int NR, int NC) { 
  std::string result = 
  "\n#define mrg31k3p_NORM_cl_T  4.656612873077392578125e-10\n\n"
  "\n\n__kernel void fisher_sim_gpu(\n"
  "   const int nrow,\n"
  "	  const int ncol,\n"
  "   __global int *nrowt, \n"
  "   __global int *ncolt, \n"
  "   const int n, \n" //ntotal
  "	  const int vsize,\n" //extra para
  "	__global double *fact,\n"
  "	__global double *results,\n" // extra para
  "	__global clrngMrg31k3pHostStream* streams"
  ") { \n"
  "	  int jwork["+ std::to_string(NR) + "];\n"  // IS THERE ENOUGH PRIVATE MEMORY TO DO THIS FOR LARGE NR?
  "	  int matrix["+ std::to_string(NR*NC) + "];\n"
  "   int i, t, u, iter, D; \n"  //original j changed to t, ii changed to u, added a D here
  "   double ans;\n"
  "   const int size = get_global_size(1)*get_global_size(0);\n"
  "   int index=get_global_id(1)*get_global_size(0) + get_global_id(0);\n"
  "   clrngMrg31k3pStream private_stream_d;\n"
  "   clrngMrg31k3pCopyOverStreamsFromGlobal(1, &private_stream_d, &streams[index]);\n"
  
  "	  int j, l, m, ia, ib, ic, jc, id, ie, ii, nll, nlm, nr_1, nc_1;\n"
  "   double x, y, dummy, sumprb;\n"
  "   bool lsm, lsp;\n"
  
  
  "   for(D = index; D < vsize; D += size) {\n"
  
  
  "      nr_1 = nrow - 1;\n"
  "      nc_1 = ncol - 1;\n"
  "      ib = 0;\n" /* -Wall */

/* Construct random matrix */
"       for (j = 0; j < nc_1; ++j)\n"
"	      jwork[j] = ncolt[j];\n"

"       jc = n;\n"
/* -----  matrix[ l, * ] ----- */
"       for (l = 0; l < nr_1; ++l) {\n"
"	    ia = nrowt[l];\n"
"	    ic = jc;\n"
"	    jc -= ia;\n"/* = n_tot - sum(nr[0:l]) */

"    for (m = 0; m < nc_1; ++m) {\n"
"    id = jwork[m];\n"
"    ie = ic;\n"
"    ic -= id;\n"
"    ib = ie - ia;\n"
"    ii = ib - id;\n"

"    if (ie == 0) {\n" /* Row [l,] is full, fill rest with zero entries */
"       for (j = m; j < nc_1; ++j)\n"
"       matrix[l + j * nrow] = 0;\n"
"       ia = 0;\n"
"       break;\n"
"    }\n" // if ie

/* Generate pseudo-random number */
"    dummy = clrngMrg31k3pNextState(&private_stream_d.current) * mrg31k3p_NORM_cl_T;\n"


"     //do { \n"// Outer Loop  COMMENTED OUT
/* Compute conditional expected value of MATRIX(L, M) */
"     nlm = (int)(ia * (id / (double) ie) + 0.5);\n"
"     x = exp(fact[ia] + fact[ib] + fact[ic] + fact[id]- fact[ie] - fact[nlm]\n"
"     - fact[id - nlm] - fact[ia - nlm] - fact[ii + nlm]);\n"
"     if (x >= dummy)\n"
"     break;\n"
//  "     if (x == 0.) \n"       /* MM: I haven't seen this anymore */
//  "     error(_(\"rcont2 [%d,%d]: exp underflow to 0; algorithm failure\"), l, m);\n"

"     sumprb = x;\n"
"     y = x;\n"
"     nll = nlm;\n"

" //do {\n"     /* Increment entry in row L, column M  COMMENTED OUT!!*/
" j = (int)((id - nlm) * (double)(ia - nlm));\n"
" lsp = (j == 0);\n"
" if (!lsp) {\n"
" ++nlm;\n"
" x = x * j / ((double) nlm * (ii + nlm));\n"
" sumprb += x;\n"
" if (sumprb >= dummy)\n"
" goto L160;\n"
" }\n" // if !lsp

#ifdef UNDEF  
" do {\n"  // do 2 
//  " R_CheckUserInterrupt();\n"

/* Decrement entry in row L, column M */
" j = (int)(nll * (double)(ii + nll));\n"
" lsm = (j == 0);\n"
" if (!lsm) {\n"
" --nll;\n"
" y = y * j / ((double) (id - nll) * (ia - nll));\n"
" sumprb += y;\n"
"  if (sumprb >= dummy) {\n"
"  nlm = nll;\n"
"  goto L160;\n"
"  }\n" // if sumprb
"  if (!lsp)\n"
"  break;\n"/* to while (!lsp) */
" }\n" // if !lsp
" } while (!lsm);\n" // do 2
#endif

"// } while (!lsp);\n" // do increment entry COMMENTED OUT!!

" dummy = sumprb * clrngMrg31k3pNextState(&private_stream_d.current) * mrg31k3p_NORM_cl_T;\n"

" //} while (1);\n"  // do outer loop  !! COMMENTED OUT!

"L160:\n"
"matrix[l + m * nrow] = nlm;\n"
"ia -= nlm;\n"
"jwork[m] -= nlm;\n"
"}\n" // for m
"matrix[l + nc_1 * nrow] = ia;\n"/* last column in row l */
"}\n" // for l

"for (m = 0; m < nc_1; ++m)\n"
"matrix[nr_1 + m * nrow] = jwork[m];\n"

"matrix[nr_1 + nc_1 * nrow] = ib - matrix[nr_1 + (nc_1-1) * nrow ];\n"


"ans = 0.;\n"
"for (t = 0; t < ncol; ++t) {\n"
"for (i = 0, u = j * nrow; i < nrow;  i++, u++)\n"
"ans -= fact[matrix[u]];\n"
"}\n"  // for t


"results[D] = ans;\n"


"}\n" // for D loop
"clrngMrg31k3pCopyOverStreamsToGlobal(1,  &streams[index], &private_stream_d);\n"
"}\n" ;
 return(result); 
}

void fisher_sim_gpu(
    viennacl::vector_base<int> &sr, 
    viennacl::vector_base<int> &sc,
    viennacl::vector_base<double> &ans, 
    Rcpp::IntegerMatrix streamsR, 
    Rcpp::IntegerVector numWorkItems,
    int ctx_id){
  
  const int nr = sr.size(), nc = sc.size();
  //  scalar<int> nScalar=-1;
  //  viennacl::linalg::sum_impl(sr, nScalar); 
  //  sum_impl(sr, nScalar);
  int n = viennacl::linalg::sum(sr);
  int vsize= ans.size();
  
  viennacl::ocl::context ctx(viennacl::ocl::get_context(ctx_id));
  
  std::string random_kernel_string = mrg31k3pTypeString<double>();
  std::string kernel_string= random_kernel_string + FisherSimkernelString(nr, nc);

#ifdef DEBUG
  Rcpp::Rcout << kernel_string << "\n\n";
#endif  
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
  
  

  


  
//  viennacl::vector_base<int> observed = viennacl::vector_base<int>(nr*nc, ctx); 

    // TO DO: put this in local memory (?)
  viennacl::vector_base<double> fact = viennacl::vector_base<double>(n+1, ctx); 
  
    /* Calculate log-factorials.  fact[i] = lgamma(i+1)
   */
  fact(0) = 0.;
  fact(1) = 0.;
  int i;
  for(i = 2; i <= n; i++) fact(i) = fact(i - 1) + log(i);
  

  viennacl::ocl::enqueue(fisher_sim(
    nr, nc, 
    sr, sc, 
    n, vsize, 
    //observed, 
    fact, 
    //jwork, 
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
