//"#include <R.h>\n"
//"#include <Rinternals.h>\n"
//"#include <Rmath.h>\n"
//"#include <R_ext/Random.h>\n"

//"#ifdef HAVE_CONFIG_H\n"
//"#include <config.h>\n"
//"#endif\n"

//"#ifdef ENABLE_NLS\n"
//"#include <libintl.h>\n"
//"#define _(String) gettext (String)\n"
//"#else\n"
//"#define _(String) (String)\n"
//"#endif\n"


//"#include <R_ext/Applic.h>\n"
//"#include <R_ext/Boolean.h>\n"
//"#include <R_ext/Error.h>\n"
//"#include <R_ext/Utils.h>\n"

#include <string>

std::string FisherSimkernelString  = 
  
  "\n\n__kernel void fisher_sim_gpu(\n"
  "const int nrow,\n"
  "	  const int ncol,\n"
  "   __global int *nrowt, \n"
  "   __global int *ncolt, \n"
  "   const int n, \n" //ntotal
  "	  int vsize,\n" //extra para
  "	__global int *matrix, \n"
 "	__global double *fact,\n"
 "	__global int *jwork, \n"
  "	__global double *results,\n" // extra para
  "	__global clrngMrg31k3pHostStream* streams"
  ") { \n"
"   int i, t, u, iter; \n"  //original j changed to t, ii changed to u
  "   double ans;\n"
  "   const int size = (get_global_size(1)*get_global_size(0));\n"
  "   int index=get_global_id(1)*get_global_size(0) + get_global_id(0);\n"
  
  /* Calculate log-factorials.  fact[i] = lgamma(i+1) */
  "   fact[0] = fact[1] = 0.;\n"
  "   for(i = 2; i <= n; i++)\n"
  "	fact[i] = fact[i - 1] + log(i);\n"
  
  "   clrngMrg31k3pStream private_stream_d;\n"
  "   clrngMrg31k3pCopyOverStreamsFromGlobal(1, &private_stream_d, &streams[index]);\n"
  
  "   for(D = index; D < vsize; D += size) {\n"
  
  "	   int j, l, m, ia, ib, ic, jc, id, ie, ii, nll, nlm, nr_1, nc_1;\n"
  "      double x, y, dummy, sumprb;\n"
  "      Rboolean lsm, lsp;\n"
  
  "      nr_1 = nrow - 1;\n"
  "      nc_1 = ncol - 1;\n"
  "      ib = 0;\n" /* -Wall */
  
  /* Construct random matrix */
  "       for (j = 0; j < nc_1; ++j)\n"
  "	      jwork[j] = ncolt[j];\n"
  
  "       jc = n;\n"
  
  "       for (l = 0; l < nr_1; ++l) {\n" /* -----  matrix[ l, * ] ----- */
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
  "    }\n"
  
  /* Generate pseudo-random number */
  "    dummy = clrngMrg31k3pNextState(&private_stream_d.current) * mrg31k3p_NORM_cl_double;\n"
  
  "     do { \n"        /* Outer Loop */
  /* Compute conditional expected value of MATRIX(L, M) */
  "     nlm = (int)(ia * (id / (double) ie) + 0.5);\n"
  "     x = exp(fact[ia] + fact[ib] + fact[ic] + fact[id]- fact[ie] - fact[nlm]\n"
  "     - fact[id - nlm] - fact[ia - nlm] - fact[ii + nlm]);\n"
  "     if (x >= dummy)\n"
  "     break;\n"
  "     if (x == 0.) \n"       /* MM: I haven't seen this anymore */
  "     error(_(\"rcont2 [%d,%d]: exp underflow to 0; algorithm failure\"), l, m);\n"
  
  "     sumprb = x;\n"
  "     y = x;\n"
  "     nll = nlm;\n"
  
  " do {\n"     /* Increment entry in row L, column M */
  " j = (int)((id - nlm) * (double)(ia - nlm));\n"
  " lsp = (j == 0);\n"
  " if (!lsp) {\n"
  " ++nlm;\n"
  " x = x * j / ((double) nlm * (ii + nlm));\n"
  " sumprb += x;\n"
  " if (sumprb >= dummy)\n"
  " goto L160;\n"
  " }\n"
  
  " do {\n"
  " R_CheckUserInterrupt();\n"
  
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
  "  }\n"
  "  if (!lsp)\n"
  "  break;\n"/* to while (!lsp) */
  " }\n"
  " } while (!lsm);\n"
  
  " } while (!lsp);\n"
  
  " dummy = sumprb * clrngMrg31k3pNextState(&private_stream_d.current) * mrg31k3p_NORM_cl_double;\n"
  
  " } while (1);\n"
  
  "L160:\n"
  "matrix[l + m * nrow] = nlm;\n"
  "ia -= nlm;\n"
  "jwork[m] -= nlm;\n"
  "}\n"
  "matrix[l + nc_1 * nrow] = ia;\n"/* last column in row l */
  "}\n"
  
  "for (m = 0; m < nc_1; ++m)\n"
  "matrix[nr_1 + m * nrow] = jwork[m];\n"
  
  "matrix[nr_1 + nc_1 * nrow] = ib - matrix[nr_1 + (nc_1-1) * nrow ];\n"
  
  
  "ans = 0.;\n"
  "for (t = 0; t < ncol; ++t) {\n"
  "for (i = 0, u = j * nrow; i < nrow;  i++, u++)\n"
  "ans -= fact[matrix[u]];\n"
  "}\n"
  
  
  "results[D] = ans;\n"
  
  
  "}\n"
  
  "clrngMrg31k3pCopyOverStreamsToGlobal(1,  &streams[index], &private_stream_d);\n"
  
  
  "}\n"
  ;