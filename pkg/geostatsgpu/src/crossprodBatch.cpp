#include "geostatsgpu.hpp"
//#define DEBUG

// C = A^T A 
// TO DO: C = A^T f(D) A


template <typename T> 
std::string crossprodBatchString(
    const int Nrow, 
    const int Ncol,
    const int Nmatrix, 
    const int NpadC, 
    const int NpadA, 
    const int NpadBetweenMatricesC,
    const int NpadBetweenMatricesA,
    const int NlocalCacheA, // numbers of rows to cache of A
    const std::vector<int> Nlocal// cache a NlocalCacheC by NlocalCacheC submatrix of C
  ) { 
    
 /*
  * global groups col by matrix
  * local items inner by row
  * 
  */
 
  std::string typeString = openclTypeString<T>();
  std::string result = "";
  
  if(typeString == "double") {
    result += "\n#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n\n";
  }
  result += 
    "\n#define Nrow " + std::to_string(Nrow) + "\n"    
    "#define Ncol " + std::to_string(Ncol) + "\n"    
    "#define Nmatrix " + std::to_string(Nmatrix) + "\n"    
    "#define NpadC " + std::to_string(NpadC) + "\n"    
    "#define NpadA " + std::to_string(NpadA) + "\n"    
    "#define NpadLocal " + std::to_string(Nlocal[1]) + "\n"    
    "#define NpadBetweenMatricesC " + std::to_string(NpadBetweenMatricesC) + "\n"    
    "#define NpadBetweenMatricesA " + std::to_string(NpadBetweenMatricesA) + "\n"    
    "#define NrowStop " + std::to_string(std::min(NlocalCacheA, Nrow)) + "\n"    
    "#define NlocalCacheA "  + std::to_string(NlocalCacheA) + "\n\n";   
    
    result +=
    "__kernel void crossprodBatch(\n"
    "	__global " + typeString+ " *C,\n"
    "	__global "+ typeString+ " *B) {\n\n"
      
    "local " + typeString + " Acache[" + 
      std::to_string(NlocalCacheA) + "];\n" 
    "local " + typeString + " Ccache[" + 
      std::to_string(Nlocal[0] * Nlocal[1]) + "];\n" +

  typeString + "Cout, Ctemp;\n"
  "event_t ev\n;"
  "int AHere, CHere;\n"
  "int Dmatrix, Drow, Dcol, DrowNpadC, DcolNpadA, DrowNpadA, Dinner;\n"
  "const int AHereInc = get_num_groups(1)*NpadBetweenMatricesA;\n"
  "const int CHereInc = get_num_groups(1)*NpadBetweenMatricesC;\n"
  "const int DlocalInc = get_local_size(0)*NpadLocal;\n"
  "const int DlocalIncDiag = DlocalInc * DrowNpadA;\n"
  "const int DcolNpadAInc = get_num_groups(0)*NpadA;\n"
  "const int DrowNpadCInc = get_local_size(1)*NpadC;\n"
  "const int DrowNpadAInc = get_local_size(1)*NpadA;\n"
  "const int doLocalSum = (get_local_id(0)==0);\n"
  "const int doFinalSum = (get_local_id(0)==0 & get_local_id(1)==0);\n"
  "const int cacheIndex = get_local_id(1)+NpadLocal*get_local_id(0);\n";
  
      result +=  "\n\n"
  "for(Dmatrix = get_group_id(1),\n"
  "    AHere = Dmatrix * NpadBetweenMatricesA,\n"
  "    CHere = Dmatrix * NpadBetweenMatricesC;\n"
  "  Dmatrix < Nmatrix;\n"
  "  Dmatrix += get_num_groups(1),\n"
  "    AHere += AHereInc,\n"
  "    CHere += CHereInc){\n";

  result +=  "\n"
  "  for(Dcol = get_group_id(0),\n"
  "      DcolNpadA = AHere + Dcol * NpadA;\n"
  "    Dcol < Ncol;\n"
  "    Dcol += get_num_groups(0),\n"
  "      DcolNpadA += DcolNpadAInc){\n";
  
  
  /*\
   * TO DO: 
   * switch Drow, Dinner loops, Nrounds cache
   */
  
  result +=  "\n"
  "    ev=async_work_group_strided_copy (Acache, &A[AHere],\n"
  "      NlocalCacheA, NpadA, 0);\n"
  "    wait_group_events (1, &ev);\n";

  // diagonal
  result += 
    "    Cout=0.0;\n"
    "    for(Dinner=cacheIndex;Dinner < NrowStop; Dinner += DlocalInc){\n"
    "      Cout += Acache[Dinner]*Acache[Dinner];\n"
    "    }\n"
    "    for(Dinner = NrowStop + cacheIndex,\n"
    "          DrowNpadA = DcolNpadA + Dinner * NpadA;\n"
    "        Dinner < Nrow;\n"
    "        Dinner += DlocalInc, DrowNpadA += DlocalIncDiag){\n"
    "      Ctemp = A[DrowNpadA];\n"
    "      Cout += Ctemp * Ctemp;\n"
    "    }\n"
    "    Ccache[cacheIndex] = Cout;\n"
    "    barrier(CLK_LOCAL_MEM_FENCE);\n"
    "    if(doFinalSum) {\n"
    "     for(Drow = 0,DrowNpadA=0; Drow < get_local_size(1);\n"
    "         Drow++,DrowNpadA += NpadLocal){\n"
    "       for(Dinner = 0;Dinner < get_local_size(0); Dinner++){\n"
    "         Cout += Ccache[DrowNpadA + Dinner];\n"
    "       }\n"
    "     }\n"
    "     C[DcolNpadC + Dcol] = Cout;\n"
    "    }\n"
    "    barrier(CLK_LOCAL_MEM_FENCE);\n";
    
    // off-diagonals
  result += 
  "    for(Drow = 1 + Dcol + get_local_id(1),\n"
  "        DrowNpadC = CHere + Drow * NpadC,\n"
  "        DrowNpadA = AHere + Drow * NpadA;\n"
  "      Drow < Ncolm1;\n"
  "      Drow += get_local_size(1),\n" 
  "        DcolNpadC += DrowNpadCInc,\n"
  "        DrowNpadA += DrowNpadAInc ){\n";
  
  result +=  "\n"
  "      Cout=0.0;\n"
  "      for(Dinner = get_local_id(0);\n"
  "        Dinner < Nstop;\n"
  "        Dinner += get_local_size(0)){\n";
  result += 
    "      Cout += A[DrowNpadA + Dinner] * Acache[Dinner]\n";
  
  "      for(Dinner < Nstop + get_local_id(0);\n"
  "        Dinner < Nrow;\n"
  "        Dinner += get_local_size(0)){\n";
  
  result += 
    "      Cout += A[DrowNpadA + Dinner] *\n"
    "        A[DcolNpadA + Dinner]\n";
    
    result += 
      "      }// Dinner\n";
    result +=       
      "      Ccache[cacheIndex] = Cout;\n"
      "      barrier(CLK_LOCAL_MEM_FENCE);\n"
      "      if(doLocalSum){"
      "for(Dinner = 1;Dinner < get_local_size(0);Dinner++){"
      "  Ccache[cacheIndex] += Ccache[cacheIndex + Dinner * NpadLocal;\n"
      "}\n"
      "C[DrowNpadC + Dcol] = Ccache[cacheIndex];\n"
      "      }\n"
      "      barrier(CLK_LOCAL_MEM_FENCE);\n";
      
    result += 
      "    }// Drow\n";
    result += 
      "  }// Dcol\n";
    result += 
      "}// Dmatrix\n";
    result += 
      "}// function";
  
  return(result);
}
