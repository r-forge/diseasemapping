//#include <CL/fisher_sim.hpp>   
#include "random_number.hpp"

//#define DEBUG

using namespace Rcpp;
using namespace viennacl;	
using namespace viennacl::linalg;


// adapted from https://github.com/wch/r-source/blob/HEAD/src/library/stats/src/rcont.c
// https://github.com/wch/r-source/blob/HEAD/src/library/stats/src/chisqsim.c
template <typename T> 
std::string FisherSimkernelString(int NR, int NC) { 
  
  std::string typeString = openclTypeString<T>();
  std::string result = mrg31k3pTypeString<T>(); 

  
  result += "\n"
  "\n#define MAXITER 10\n"
  "#define nrow " + std::to_string(NR) + "\n"
  "#define ncol " + std::to_string(NC) + "\n"
  "#define nr_1  " + std::to_string(NR - 1) + "\n"
  "#define nc_1 " + std::to_string(NC - 1) + "\n"
  "#define nc_1nrow " + std::to_string((NC - 1) * NR) + "\n"
  "#define nrowncol " + std::to_string(NR*NC) + "\n"
  "\n__kernel void fisher_sim_gpu(\n"

  //  "   const int nrow,\n"
//  "	  const int ncol,\n"
  "   __global int *nrowt, \n"
  "   __global int *ncolt, \n"
  "   const int n, \n" //ntotal
  "	  const int vsize,\n" //extra para
  "	__global  " + typeString + "  *fact,\n"
  "	__global  " + typeString + "  *results,\n" // fisher log p
#ifdef DEBUGEXTRA  
  "	__global  " + typeString + "  *extraR,\n" // random numbers
  "	__global  " + typeString + "  *extraX,\n" // x statistics
#endif  
  "	__global clrngMrg31k3pHostStream* streams"
  ") { \n\n"
  "	  local int jwork[ncol];\n"
  "	  local int matrix[nrowncol];\n"
  "   int i, t, u, iter, D; \n"  //original j changed to t, ii changed to u, added a D here
  "   " + typeString + "  ans;\n"
  "   const int globalSize = get_global_size(1)*get_global_size(0);\n"
  "   const int index=get_global_id(1)*get_global_size(0) + get_global_id(0);\n"
  
  "	  local int jjj, l, m, ia, ib, ic, jc, id, ie, ii, nll, nlm;\n"
  "   local " + typeString + " x, y, dummy, sumprb;\n"
  "   bool lsm, lsp;\n\n"
  "   int Diter, Diter2, Diter3, goTo160;\n"

  "   clrngMrg31k3pStream private_stream_d;\n"
  "   clrngMrg31k3pCopyOverStreamsFromGlobal(1, &private_stream_d, &streams[index]);\n\n"
  
  
"   for(D = index; D < vsize; D += globalSize) {\n"

  
"      ib = 0;\n" /* -Wall */
  
/* Construct random matrix */
#ifdef UNDEF
"      for (jjj = 0; jjj < nc_1; ++jjj){\n"
"	       jwork[jjj] = ncolt[jjj];\n"
//"  extraR[jj + D*nrowncol*size+Diter*nrowncol*size*vsize]= jwork[jjj];"
"      }\n"
#endif
"async_work_group_copy(jwork, ncolt, nc_1, 0);\n"


"      jc = n;\n"
/* -----  matrix[ l, * ] ----- */
"      for (l = 0; l < nr_1; ++l) {\n"
"	         ia = nrowt[l];\n"
"	         ic = jc;\n"
"	         jc -= ia;\n"/* = n_tot - sum(nr[0:l]) */

"          for (m = 0; m < nc_1; ++m) {\n"
"            id = jwork[m];\n"
"            ie = ic;\n"
"            ic -= id;\n"
"            ib = ie - ia;\n"
"            ii = ib - id;\n"

"            if (ie == 0) {\n" /* Row [l,] is full, fill rest with zero entries */
"              for (jjj = m; jjj < nc_1; ++jjj)\n"
"                matrix[l + jjj * nrow] = 0;\n"
"              ia = 0;\n"
"              m = nc_1;// induce break of m loop\n"  
#ifdef DEBUGEXTRA
"     extraX[l + m*nrow + D*nrowncol+Diter*nrowncol*vsize]= -1.0;\n"
#endif  
// used to have "             break;\n" 
"            } else {\n" // if ie not zero

/* Generate pseudo-random number */
"            dummy = mrg31k3p_NORM_cl * clrngMrg31k3pNextState(&private_stream_d.current);\n"

// LOOP 1
"goTo160 = 0;\n"
"Diter = 0;\n"
"do { \n"
"  Diter ++;\n"
/* Compute conditional expected value of MATRIX(L, M) */
"     nlm = (int)(ia * (id / (" + typeString + ") ie) + 0.5);\n"
"     x = exp(fact[ia] + fact[ib] + fact[ic] + fact[id]- fact[ie] - fact[nlm]\n"
"         - fact[id - nlm] - fact[ia - nlm] - fact[ii + nlm]);\n"
// l is row, m column
#ifdef DEBUGEXTRA
"     extraX[l + m*nrow + D*nrowncol+Diter*nrowncol*vsize]= x;\n"
#endif  

"     if (x >= dummy) {\n"
"        break;\n"
"     }\n" // break loop 1

"     sumprb = x;\n"
"     y = x;\n"
"     nll = nlm;\n"

 // LOOP 2
" Diter2 = 0;\n"
" do {\n"
" Diter2 ++;\n"
"// j = (int)((id - nlm) * (" + typeString + ")(ia - nlm));\n"
" jjj=(int)((id - nlm) * (" + typeString + ")(ia - nlm));\n"
" lsp = (jjj == 0);\n"
" if (!lsp) {\n"
"   ++nlm;\n"
"   x = x * jjj / ((" + typeString + ") nlm * (ii + nlm));\n"
"   sumprb += x;\n"
"   if (sumprb >= dummy) {\n"
//"     goto L160;\n"
"     goTo160 = 1;\n"
"   }\n"
" }\n" // if !lsp
" if(goTo160) {"
"    break;\n" // break loop 2
" }\n"


// LOOP 3
" Diter3 = 0;\n"
" do {\n"  
" Diter3 ++;\n"
/* Decrement entry in row L, column M */
" jjj = (int)(nll * (" + typeString + ")(ii + nll));\n"  
" lsm = (jjj == 0);\n"

" if (!lsm) {\n"
" --nll;\n"
" y = y * jjj / ((" + typeString + ") (id - nll) * (ia - nll));\n"
" sumprb += y;\n"
"  if (sumprb >= dummy) {\n"
"  nlm = nll;\n"
//"  goto L160;\n"
"    goTo160 = 1;\n"
"  }\n" // if sumprb
"  if (!lsp)\n" 
"    break;\n" // to while (!lsp) loop 3
"  }\n" // if !lsm
" } while (!lsm & (Diter3 < MAXITER) & (!goTo160) );\n" // do loop 3 !!!
// END LOOP 3



" } while (!lsp & (Diter2 < MAXITER) & (!goTo160));\n" // do loop 2 !!! 
// END LOOP 2

//" dummy = sumprb * mrg31k3p_NORM_cl * clrngMrg31k3pNextState(&private_stream_d.current);\n"
//" extraR[D + Diter*vsize] = dummy;\n"    


" } while (!goTo160 &  (Diter < MAXITER) );\n"  // do loop 1
// END LOOP 1


// still in M loop
//"L160:\n"
"        matrix[l + m * nrow] = nlm;\n"
"        ia -= nlm;\n"
"        jwork[m] -= nlm;\n\n"

"        }// IF ie not zero\n" 

"     }// M LOOP\n" 
//"     // originally was       matrix[l + nc_1 * nrow] = ia;\n"/* last column in row l */
"     matrix[l + nc_1nrow] = ia;\n"/* last column in row l */
"  }// L LOOP\n" 
    /* Compute entries in last row of MATRIX */

"  for (m = 0; m < nc_1; ++m) {\n"
"    matrix[nr_1 + m * nrow] = jwork[m];\n"
"  }\n"

"  matrix[nr_1 + nc_1nrow] = ib - matrix[nr_1 + (nc_1-1) * nrow ];\n"
// end of rcont2
// now on fisher_sim
#ifdef DEBUGEXTRA
"      for (l = 0; l < nrow; ++l) {\n"
"         for (m = 0; m < ncol; ++m) {\n"
"  extraR[l + m*nrow + D*nrowncol]= matrix[nr_1 + m * nrow];\n"
"}}\n"
#endif  

"  ans = 0.;\n"
"  for (m = 0; m < ncol; ++m) {\n"
"    for (l = 0,u=m*nrow; l < nrow;  l++, u++) {\n"
"      ans -= fact[matrix[u]];\n"
#ifdef DEBUGEXTRA
"  extraR[u + D*nrowncol + 1*nrowncol*vsize]= ans;\n"
"  extraR[u + D*nrowncol + 2*nrowncol*vsize]= fact[matrix[u]];\n"
"  extraR[u + D*nrowncol + 3*nrowncol*vsize]= matrix[u];\n"
"  extraR[u + D*nrowncol + 4*nrowncol*vsize]= u;\n"
#endif  
"    }\n"  // for l
"  }\n\n"  // for m



"  results[D] = ans;\n"    


"}\n" // for D loop
"clrngMrg31k3pCopyOverStreamsToGlobal(1,  &streams[index], &private_stream_d);\n"
"}\n" ;
 return(result); 
}

void fisher_sim_gpu(
    viennacl::vector_base<int> &sr, 
    viennacl::vector_base<int> &sc,
    viennacl::vector_base<double> &ans, 
//    viennacl::vector_base<double> &extraR, 
//    viennacl::vector_base<double> &extraX, 
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
  
  std::string kernel_string= FisherSimkernelString<double>(nr, nc);

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
  fisher_sim.local_work_size(0, 1L);
  fisher_sim.local_work_size(1, 1L);
  

  
//  viennacl::vector_base<int> observed = viennacl::vector_base<int>(nr*nc, ctx); 

    // TO DO: put this in local memory (?)
  viennacl::vector_base<double> fact = viennacl::vector_base<double>(n+1, ctx); 
  
    /* Calculate log-factorials.  fact[i] = lgamma(i+1)
   */
  fact(0) = 0.;
  fact(1) = 0.;
  int i;
  for(i = 2; i <= n; i++) {
    fact(i) = fact(i - 1) + log(i);
  }


    viennacl::ocl::enqueue(fisher_sim(
//    nr, nc, 
    sr, sc, 
    n, vsize, 
    //observed, 
    fact, 
    //jwork, 
    ans,
//    extraR, extraX,
    bufIn
    ) ); //streams, out, vector_size

  // copy streams back to cpu
  viennacl::backend::memory_read(bufIn.handle(), 0, streamBufferSize, streams);
  //#endif 
  // then transfer to R object, //return streams to R 
  convertclRngMat(streams, streamsR);
}


//    Rcpp::S4  extraRR,
//    Rcpp::S4  extraXR,

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
#ifdef DEBUGEXTRA
  std::shared_ptr<viennacl::vector_base<double> > extraR = getVCLVecptr<double>(extraRR.slot("address"), 1, ctx_id);
  std::shared_ptr<viennacl::vector_base<double> > extraX = getVCLVecptr<double>(extraXR.slot("address"), 1, ctx_id);
#endif
  
  
  fisher_sim_gpu(*sr, *sc, *ans, 
//    *extraR, *extraX, 
    streamsR, max_global_size, ctx_id);
  
  return (ansR);
}
