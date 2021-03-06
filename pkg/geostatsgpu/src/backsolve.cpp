#include "geostatsgpu.hpp"
#define DEBUG

   
/*
  solve L a = b for a
  solve A C = B for C
 

 
 - Dgroup[0] is matrix Dlocal[0] is row of A and C
 - Dgroup[1] is a column of B, Dlocal[1] is Dinner
 - each Dgroup[1] does NcolsGroup columns of B
 - NcolsPerGroup = ceil( Ncol/Ngroup[1])
 - local cache Nlocal[1] submatrices of size Nlocal[0] by NcolsPerGroup of C
 - don't cache A ?
 
 - loop Dmatrix = Dgroup[0]
 - loop Drow = Dlocal[0]
 - loop Dinner = Dlocal[1]
 - cache A ?
 - loop DcolB = Dgroup[1]
 - compute part of Linv b
 - end Dinner, DcolB, 
 - sum Linv b over Dlocal[1] in cache
 - the Dlocal[1]=0's fill in a's
 
  

 */

//  solve A C = B for C

template <typename T> 
std::string backsolveBatchString(
    const int sameB,
    const int diagIsOne,
    const int Nrow, 
    const int Ncol,
    const int Nmatrix, 
    const int NpadC, 
    const int NpadA, 
    const int NpadB,
    const int NpadBetweenMatricesC,
    const int NpadBetweenMatricesA,
    const int NpadBetweenMatricesB,
    const int NrowsToCache, const int NcolsPerGroup,
    const int NlocalCacheC,  // NrowsToCache by NcolsPerGroup
    const int NlocalCacheSum, // Nlocal(0) * Nlocal(1) * NcolsPerGroup 
    const int NpadBetweenMatricesSum // Nlocal(0) * Nlocal(1)
) {
  
  std::string typeString = openclTypeString<T>();
  std::string result = "";
  
  if(typeString == "double") {
    result += "\n#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n\n";
  }
  result = result + 
    "\n#define Nrow " + std::to_string(Nrow) + "\n"    
    "#define Ncol " + std::to_string(Ncol) + "\n"    
    "#define Nmatrix " + std::to_string(Nmatrix) + "\n"    
    "#define NpadC " + std::to_string(NpadC) + "\n"
    "#define NpadA " + std::to_string(NpadA) + "\n"    
    "#define NpadB " + std::to_string(NpadB) + "\n"
    "#define NpadCcache " + std::to_string(NcolsPerGroup) + "\n"
    "#define NpadBetweenMatricesC " + std::to_string(NpadBetweenMatricesC) + "\n"    
    "#define NpadBetweenMatricesA " + std::to_string(NpadBetweenMatricesA) + "\n"    
    "#define NpadBetweenMatricesB " + std::to_string(NpadBetweenMatricesB) + "\n"    
    "#define NrowsToCache " + std::to_string(NrowsToCache) + "\n"
    "#define NrowStop " + std::to_string(std::min(Nrow, 
                                 NrowsToCache)) + "\n"
    "#define NlocalCacheC " + std::to_string(NlocalCacheC) + "\n"    
    "#define NlocalCacheSum " + std::to_string(NlocalCacheSum) + "\n"
    "#define NpadBetweenMatricesSum " + std::to_string(NpadBetweenMatricesSum) + "\n\n";
  
  result += "__kernel void backsolveBatch(\n"
  "	__global " + typeString+ " *C,\n"
  "	__global "+ typeString+ " *A,\n"
  "	__global "+ typeString+ " *B) {\n\n"
  "local int AHere, BHere, CHere, DrowZero;\n"
  "int AHereRow, BHereRow, CHereRow, CcacheHereRow;\n"
  "local "+ typeString+ " Ccache[NlocalCacheC];\n" 
  "local "+ typeString+ " cacheSum[NlocalCacheSum];\n"
  + typeString + " Acache;\n" 
  "int Dmatrix, Drow, Dcol, Dinner, DinnerC, DcolCache, Dsum;\n"
  "local int DCrowInc, DArowInc, DBrowInc, DinnerCinc, DCcacheRowInc;\n"
  "const int DlocalCache = get_local_id(0) * get_local_size(1) + get_local_id(1);\n"
  "const int localIsFirstCol = (get_local_id(1) == 0);\n"
  "const int localIsFirstItem = (get_local_id(0) == 0) & (get_local_id(1) == 0);\n";

    
  result += 
    "if (localIsFirstItem) {\n"
    "  DCrowInc = get_local_size(0) * NpadC;\n"
    "  DCcacheRowInc = get_local_size(0) * NpadCcache;\n"
    "  DArowInc = get_local_size(0) * NpadA;\n"
    "  DBrowInc = get_local_size(0) * NpadB;\n"
//    "  DinnerCinc = get_local_size(1) * NpadC;\n"
    "  DinnerCinc = get_local_size(1) * NpadCcache;\n"
    "}\n";
    
// loop through matrix
  result +=  "\n"
    "for(Dmatrix = get_group_id(0); Dmatrix < Nmatrix; Dmatrix += get_num_groups(1)){\n"
    "  barrier(CLK_LOCAL_MEM_FENCE);\n"
    "  if(localIsFirstItem) {\n"
    "    AHere = Dmatrix*NpadBetweenMatricesA;\n"
    "    CHere = Dmatrix*NpadBetweenMatricesC;\n"
    "    BHere = Dmatrix*NpadBetweenMatricesB;\n"
    "  };\n"
    "  barrier(CLK_LOCAL_MEM_FENCE);\n";

  // loop through rows of A
  result +=
    "  for(Drow = get_local_id(0),\n"
    "      BHereRow = BHere + get_local_id(0) * NpadB,\n"
    "      CHereRow = CHere + get_local_id(0) * NpadC,\n"
    "      CcacheHereRow = get_local_id(0) * NpadCcache,\n"
    "      AHereRow = AHere + get_local_id(0) * NpadA;\n"
    "      Drow < NrowStop;\n" 
    "      Drow += get_local_size(0), \n"
    "      BHereRow += DBrowInc,\n"
    "      AHereRow += DArowInc, CHereRow += DCrowInc,\n"
    "      CcacheHereRow += DCcacheRowInc){\n"

    "    barrier(CLK_LOCAL_MEM_FENCE);\n"
    "    if(localIsFirstItem) {\n"
    "      DrowZero = Drow;\n"
    "    };\n"
    "    barrier(CLK_LOCAL_MEM_FENCE);\n";

  result += 
    "    // initialize cacheSum\n"
    "    for(Dcol = get_group_id(1), DcolCache=DlocalCache;\n"
    "      Dcol < Ncol;\n"
    "      Dcol += get_num_groups(1), DcolCache+=NpadBetweenMatricesSum){\n";

//    result += "cacheSum[DcolCache] = 0.0;\n";
    result += "cacheSum[DcolCache] = 0.0;\n";
      result += 
    "    }//for Dcol\n"
    "    barrier(CLK_LOCAL_MEM_FENCE);\n";
      
      result += 
        "    for(Dinner = get_local_id(1),\n"
        "      DinnerC = get_local_id(1)*NpadCcache;\n"
        "      Dinner < DrowZero;\n"
        "      Dinner += get_local_size(1), DinnerC += DinnerCinc){\n";
      
      result += 
        "      Acache = A[AHereRow + Dinner];\n"; 
//"      Acache = A[AHere + Drow*NpadA + Dinner];\n"; 
      
  result += 
    "      // loop through columns of B and C\n"
    "      for(Dcol = get_group_id(1),DcolCache=0, Dsum=DlocalCache;\n"
    "        Dcol < Ncol;\n"
    "        Dcol += get_num_groups(1), DcolCache++, Dsum+=NpadBetweenMatricesSum){\n";


  result += 
    "        cacheSum[Dsum] += Acache * Ccache[DinnerC];\n";
  //    "        cacheSum[DlocalCache + DcolCache * NpadBetweenMatricesSum] +=\n"
  //  "           Acache * Ccache[Dinner * NpadCcache + DcolCache];\n";
  
  
  result += 
    "      }//for Dcol\n";


    result += 
    "    }//for Dinner\n"
    "    barrier(CLK_LOCAL_MEM_FENCE);\n";
    
    
    result += 
      "    // loop columns again\n"
      "    for(Dcol = get_group_id(1),DcolCache=0;\n"
      "      Dcol < Ncol;\n"
      "      Dcol += get_num_groups(1), DcolCache++){\n";
    
    result += "\n"
    "      if(localIsFirstCol){\n";
    
    result += 
      "        for(Dinner = 1, DinnerC = DlocalCache+1;\n"
      "          Dinner < get_local_size(1);\n"
      "          Dinner++,DinnerC++){\n"
      
      "            cacheSum[NpadBetweenMatricesSum*DcolCache + DlocalCache] +=\n"
      "              cacheSum[NpadBetweenMatricesSum*DcolCache + DinnerC];\n"
      
      "        }//Dinner\n"
      "        barrier(CLK_LOCAL_MEM_FENCE);\n";
    
  result +=
    "       // last bit of the triangle\n"    
    "       for(Dinner = 0,DinnerC = 0;\n"
    "         Dinner < get_local_size(0);\n"
    "         Dinner++,DinnerC += get_local_size(1)){\n";

  result +=
    "        // create C in cache and copy to C\n"
    "        if(get_local_id(0) == Dinner){\n"
    
    "          cacheSum[NpadBetweenMatricesSum*DcolCache + DinnerC] = B[BHereRow + Dcol] -\n"
    "               cacheSum[NpadBetweenMatricesSum*DcolCache + DinnerC];\n"
  
    "          C[CHereRow + Dcol] = cacheSum[NpadBetweenMatricesSum*DcolCache + DinnerC];\n"

"          Ccache[CcacheHereRow + DcolCache] =\n"
//"          Ccache[Drow*NpadCcache + DcolCache] =\n"
"            cacheSum[NpadBetweenMatricesSum*DcolCache + DinnerC];\n"
    "        } //if(get_local_id(0) == Dinner)\n"
    "        barrier(CLK_LOCAL_MEM_FENCE);\n"; 
  

    result +=
    "        // update A[Drow, ] * B[, Dcol] \n"
    "        if(get_local_id(0) > Dinner){\n"
    "          cacheSum[NpadBetweenMatricesSum*DcolCache + DlocalCache] +=\n"
"            A[AHereRow + DrowZero + Dinner] * cacheSum[NpadBetweenMatricesSum*DcolCache + DinnerC];\n"
"        }//if(get_local_id(0) > Dinner)\n"
    "        barrier(CLK_LOCAL_MEM_FENCE);\n"; 

    
      result +=
        "        }//Dinner\n";


      result +=
    "      }//localIsFirstCol\n"
    "      barrier(CLK_LOCAL_MEM_FENCE);\n";  
    
  result += 
    "    }//for Dcol\n";
  result += 
    "  }//for Drow\n";

  result += 
    "  //Now rows that aren't all cached\n";
  
  result +=
    "  for(Drow = NrowStop + get_local_id(0),\n"
    "      BHereRow = BHere + Drow * NpadB,\n"
    "      CHereRow = CHere + Drow * NpadC,\n"
    "      CcacheHereRow = Drow * NpadCcache,\n"
    "      AHereRow = AHere + Drow * NpadA;\n"
    "      Drow < Nrow;\n" 
    "      Drow += get_local_size(0), \n"
    "      BHereRow += DBrowInc,\n"
    "      AHereRow += DArowInc, CHereRow += DCrowInc,\n"
    "      CcacheHereRow += DCcacheRowInc){\n"
    
    "    barrier(CLK_LOCAL_MEM_FENCE);\n"
    "    if(localIsFirstItem) {\n"
    "      DrowZero = Drow;\n"
    "    };\n"
    "    barrier(CLK_LOCAL_MEM_FENCE);\n";
  
  result += 
    "    // initialize cacheSum\n"
    "    for(Dcol = get_group_id(1), DcolCache=DlocalCache;\n"
    "      Dcol < Ncol;\n"
    "      Dcol += get_num_groups(1), DcolCache+=NpadBetweenMatricesSum){\n";
  
  //    result += "cacheSum[DcolCache] = 0.0;\n";
  result += "cacheSum[DcolCache] = 0.0;\n";
  result += 
    "    }//for Dcol\n"
    "    barrier(CLK_LOCAL_MEM_FENCE);\n";
  
  result += 
    "    for(Dinner = get_local_id(1),\n"
    "      DinnerC = get_local_id(1)*NpadCcache;\n"
    "      Dinner < NrowStop;\n"
    "      Dinner += get_local_size(1), DinnerC += DinnerCinc){\n";
  
  result += 
    "      Acache = A[AHereRow + Dinner];\n"; 
  //"      Acache = A[AHere + Drow*NpadA + Dinner];\n"; 
  
  result += 
    "      // loop through columns of B and C\n"
    "      for(Dcol = get_group_id(1),DcolCache=0, Dsum=DlocalCache;\n"
    "        Dcol < Ncol;\n"
    "        Dcol += get_num_groups(1), DcolCache++, Dsum+=NpadBetweenMatricesSum){\n";
  
  
  result += 
    "        cacheSum[Dsum] += Acache * Ccache[DinnerC];\n";
  //    "        cacheSum[DlocalCache + DcolCache * NpadBetweenMatricesSum] +=\n"
  //  "           Acache * Ccache[Dinner * NpadCcache + DcolCache];\n";
  
  
  result += 
    "      }//for Dcol\n";
  
  
  result += 
    "    }//for Dinner\n"
    "    barrier(CLK_LOCAL_MEM_FENCE);\n";
  
  result += 
    "    for(Dinner = NrowStop + get_local_id(1),\n"
    "      DinnerC = Dinner*NpadCcache;\n"
    "      Dinner < Nrow;\n"
    "      Dinner += get_local_size(1), DinnerC += DinnerCinc){\n";
  
  result += 
    "      Acache = A[AHereRow + Dinner];\n"; 
  //"      Acache = A[AHere + Drow*NpadA + Dinner];\n"; 
  
  result += 
    "      // loop through columns of B and C\n"
    "      for(Dcol = get_group_id(1),DcolCache=0, Dsum=DlocalCache;\n"
    "        Dcol < Ncol;\n"
    "        Dcol += get_num_groups(1), DcolCache++, Dsum+=NpadBetweenMatricesSum){\n";
  
  
  result += 
    "        cacheSum[Dsum] += Acache * C[CHere + Dcol + NpadC * Dinner];\n";
  //    "        cacheSum[DlocalCache + DcolCache * NpadBetweenMatricesSum] +=\n"
  //  "           Acache * Ccache[Dinner * NpadCcache + DcolCache];\n";
  
  
  result += 
    "      }//for Dcol\n";
  
  
  result += 
    "    }//for Dinner\n"
    "    barrier(CLK_LOCAL_MEM_FENCE);\n";
  
  result += 
    "    // loop columns again\n"
    "    for(Dcol = get_group_id(1),DcolCache=0;\n"
    "      Dcol < Ncol;\n"
    "      Dcol += get_num_groups(1), DcolCache++){\n";
  
  result += "\n"
  "      if(localIsFirstCol){\n";
  
  result += 
    "        for(Dinner = 1, DinnerC = DlocalCache+1;\n"
    "          Dinner < get_local_size(1);\n"
    "          Dinner++,DinnerC++){\n"
    
    "            cacheSum[NpadBetweenMatricesSum*DcolCache + DlocalCache] +=\n"
    "              cacheSum[NpadBetweenMatricesSum*DcolCache + DinnerC];\n"
    
    "        }//Dinner\n"
    "        barrier(CLK_LOCAL_MEM_FENCE);\n";
  
  result +=
    "       // last bit of the triangle\n"    
    "       for(Dinner = 0,DinnerC = 0;\n"
    "         Dinner < get_local_size(0);\n"
    "         Dinner++,DinnerC += get_local_size(1)){\n";
  
  result +=
    "        // create C in cache and copy to C\n"
    "        if(get_local_id(0) == Dinner){\n"
    
    "          cacheSum[NpadBetweenMatricesSum*DcolCache + DinnerC] = B[BHereRow + Dcol] -\n"
    "               cacheSum[NpadBetweenMatricesSum*DcolCache + DinnerC];\n"
    
    "          C[CHereRow + Dcol] = cacheSum[NpadBetweenMatricesSum*DcolCache + DinnerC];\n"
    
    "        } //if(get_local_id(0) == Dinner)\n"
    "        barrier(CLK_LOCAL_MEM_FENCE);\n"; 
  
  
  result +=
    "        // update A[Drow, ] * B[, Dcol] \n"
    "        if(get_local_id(0) > Dinner){\n"
    "          cacheSum[NpadBetweenMatricesSum*DcolCache + DlocalCache] +=\n"
    "            A[AHereRow + DrowZero + Dinner] * cacheSum[NpadBetweenMatricesSum*DcolCache + DinnerC];\n"
    "        }//if(get_local_id(0) > Dinner)\n"
    "        barrier(CLK_LOCAL_MEM_FENCE);\n"; 
  
  
  result +=
    "        }//Dinner\n";
  
  
  result +=
    "      }//localIsFirstCol\n"
    "      barrier(CLK_LOCAL_MEM_FENCE);\n";  
  
  result += 
    "    }//for Dcol\n";
  result += 
    "  }//for Drow\n";
  
  result += 
    "}\n" // Dmatrix
    "}";
  /*
   * TO DO: do remaining rows
   * multiply A with cached C
   *   copy C over
   * loop through cached blocks of C
   * do triangle bit
   *    
   */
  return(result);
}



template <typename T> 
void backsolveBatch(
    viennacl::matrix<T> &C,
    viennacl::matrix<T> &A,
    viennacl::matrix<T> &B,
    const int diagIsOne,
    std::vector<int> Nglobal,
    std::vector<int> Nlocal,
    const int NlocalCache,
    const int ctx_id) {
  
  
  const int Nrow = A.size2(), Ncol = B.size2();
  const int Nmatrix = C.size1()/Nrow;
  const int Ngroups1 = static_cast<T>(Nglobal[1]) / 
    static_cast<T>(Nlocal[1]);
  const int NcolsPerGroup = std::ceil( static_cast<T>(Ncol) / 
                                 static_cast<T>(Ngroups1));
  const int NrowsToCache = std::floor(
    static_cast<T>(NlocalCache) /static_cast<T>(NcolsPerGroup));
  
  
  // the context
  viennacl::ocl::context ctx(viennacl::ocl::get_context(ctx_id));
  
  cl_device_type type_check = ctx.current_device().type();
  
  std::string clString =  backsolveBatchString<T>(  
    Nrow == B.size1(),
    diagIsOne,
    Nrow, 
    Ncol, // ncol
    Nmatrix,
    C.internal_size2(), 
    A.internal_size2(), 
    B.internal_size2(),
    C.internal_size2()*Nrow,//NpadBetweenMatricesC,
    A.internal_size2()*Nrow,//NpadBetweenMatricesA,
    B.internal_size2()*Nrow,//NpadBetweenMatricesB,
    NrowsToCache, NcolsPerGroup, 
    NcolsPerGroup * NrowsToCache,//    NlocalCacheC, 
    Nlocal[0] * Nlocal[1] * NcolsPerGroup, //    NlocalCacheSum 
    Nlocal[0] * Nlocal[1] //NpadBetweenMatricesSum
    );  //NcolsInCache
  
#ifdef DEBUG
  
  Rcpp::Rcout << clString << "\n\n";
  
#endif  
  
  
  
  viennacl::ocl::program & my_prog = ctx.add_program(
    clString, "my_kernel");
  
  viennacl::ocl::kernel & backsolveKernel = my_prog.get_kernel("backsolveBatch");
  
  backsolveKernel.global_work_size(0, Nglobal[0]);
  backsolveKernel.global_work_size(1, Nglobal[1]);

  backsolveKernel.local_work_size(0, Nlocal[0]);
  backsolveKernel.local_work_size(1, Nlocal[1]);


  viennacl::ocl::enqueue(backsolveKernel(
      C, A, B));
  
}



template <typename T> 
SEXP backsolveBatchTyped(
    Rcpp::S4 CR,
    Rcpp::S4 AR,
    Rcpp::S4 BR,
    const int diagIsOne,
    Rcpp::IntegerVector NglobalR,
    Rcpp::IntegerVector NlocalR, 
    const int NlocalCache) {
  
  std::vector<int> Nglobal = Rcpp::as<std::vector<int> >(NglobalR);
  std::vector<int> Nlocal = Rcpp::as<std::vector<int> >(NlocalR);
  
  const int ctx_id = INTEGER(CR.slot(".context_index"))[0]-1;
  const bool BisVCL=1;
  
  
  
  std::shared_ptr<viennacl::matrix<T> > 
    AG = getVCLptr<T>(AR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > 
    BG = getVCLptr<T>(BR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > 
    CG = getVCLptr<T>(CR.slot("address"), BisVCL, ctx_id);

  backsolveBatch<T>(*CG, *AG, *BG, 
                    diagIsOne, 
                    Nglobal, Nlocal, NlocalCache, ctx_id);
  
  return Rcpp::wrap(0L);
  
}



// [[Rcpp::export]]
SEXP backsolveBatchBackend(
    Rcpp::S4 C,
    Rcpp::S4 A,
    Rcpp::S4 B,
    const int diagIsOne,
    Rcpp::IntegerVector Nglobal,
    Rcpp::IntegerVector Nlocal, 
    const int NlocalCache) {
  
  SEXP result;
  
  Rcpp::traits::input_parameter< std::string >::type classVarR(RCPP_GET_CLASS(C));
  std::string precision_type = (std::string) classVarR;
  
  if(precision_type == "fvclMatrix") {
    result = backsolveBatchTyped<float>(
      C, A, B, diagIsOne, Nglobal, Nlocal, NlocalCache);
  } else if (precision_type == "dvclMatrix") {
    result = backsolveBatchTyped<double>(
      C, A, B, diagIsOne, Nglobal, Nlocal, NlocalCache);
  } else {
    result = Rcpp::wrap(1L);
  }
  return(result);
  
}

