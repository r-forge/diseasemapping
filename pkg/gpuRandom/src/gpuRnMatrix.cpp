#include "gpuRandom.hpp"

//#define DEBUG
//#define DEBUGKERNEL



using namespace Rcpp;
using namespace viennacl;	
using namespace viennacl::linalg;




template <typename T> 
std::string mrg31k3pMatrixString(
    const int Nrow, 
    const int Ncol,
    const int NpadCol,
    const int NpadStreams,
    const std::string random_type) { 
  
  std::string typeString = openclTypeString<T>();
  std::string result = "";
  
  if(typeString == "double") {
    
    result += "\n#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n"
    "#define TWOPI 6.283185307179586 \n" 
    "#define mrg31k3p_NORM_cl 4.656612873077392578125e-10\n"
    "//TWOPI * mrg31k3p_NORM_cl\n"
    "#define TWOPI_mrg31k3p_NORM_cl 2.925836158534319248049e-09\n\n";
    
  } else if(typeString == "float") {
    result += "\n#define TWOPI 6.2831853\n" 
    "\n#define mrg31k3p_NORM_cl 4.6566126e-10\n\n"
    "//TWOPI * mrg31k3p_NORM_cl\n"
    "#define TWOPI_mrg31k3p_NORM_cl 2.9258361e-09\n\n";
  } else if (typeString == "int"){
    result += 
      "\n#define mrg31k3p_NORM_cl 1L\n";    
  } 
  // else if (random_type == "exponential"){
  //   result +=
  //  "#define mrg31k3p_NORM_cl 4.656612875245796923096e-10\n";
  // }
  
  result += 
    "\n#define mrg31k3p_M1 2147483647\n"             /* 2^31 - 1 */
  "#define mrg31k3p_M2 2147462579\n"             /* 2^31 - 21069 */
  
  "#define mrg31k3p_MASK12 511  \n"              /* 2^9 - 1 */
  "#define mrg31k3p_MASK13 16777215  \n"         /* 2^24 - 1 */
  "#define mrg31k3p_MASK2 65535     \n"          /* 2^16 - 1 */
  "#define mrg31k3p_MULT2 21069\n";
  
  result +=
    "\n#define Nrow " + std::to_string(Nrow) + "\n"    
    "#define Ncol " + std::to_string(Ncol) + "\n"
    "#define NpadStreams " + std::to_string(NpadStreams) + "\n"
    "#define NpadCol " + std::to_string(NpadCol) + "\n";    
  
  result += mrg31k3pString();
  
  result += 
    "\n\n__kernel void mrg31k3pMatrix(\n"
    "  __global int* streams,\n" 
    "  __global " + typeString + "* out){\n\n";  
  
  
  result += 
    "const int index = get_global_id(0)*get_global_size(1) + get_global_id(1);\n";
  
  result += "int Drow, Dcol, DrowStart, Dentry;\n";
  
  result += "const int DrowInc = 2*get_global_size(0), DrowStartInc = DrowInc * NpadCol;\n";
  
  result += "uint temp;\n";
  result += "uint g1[3], g2[3];\n"
  "const int startvalue=index * NpadStreams;\n";
  
  
  if(random_type == "normal"){  
    result += 
      "local " + typeString + "  part[2], cosPart1, sinPart1;\n";// local size must be 1,2
  }
  
  
  /*result +=  " for(Drow = 0, DrowStart = startvalue, Dcol = DrowStart + 3;\n"
   "   Drow < 3; Drow++, DrowStart++, Dcol++){\n"
   "   g1[Drow] = streams[DrowStart];\n"
   "   g2[Drow] = streams[Dcol];\n"
   " }\n";    */
  
  result += "streamsToPrivate(streams,g1,g2,startvalue);\n";
  
  // Drow, Dcol = 2* global(0), group(1)
  // local work item (0,0) does out[Drow, Dcol]
  // local work item (0,1) does out[Drow+1, Dcol]
  result += 
    "for(Drow=2*get_global_id(0), DrowStart = (Drow + get_local_id(1))* NpadCol;\n" 
    "    Drow < Nrow; Drow += DrowInc, DrowStart += DrowStartInc) {\n";
  
  result += 
    "  for(Dcol=get_group_id(1), Dentry = DrowStart + Dcol;\n" 
    "      Dcol < Ncol;\n" 
    "      Dcol += get_num_groups(1), Dentry += get_num_groups(1)) {\n"
    
    "    temp = clrngMrg31k3pNextState(g1, g2);\n";
  
  /* result += 
   
   // first component
   "	y1 = ((g1[1] & mrg31k3p_MASK12) << 22) + (g1[1] >> 9)\n"
   "		+ ((g1[2] & mrg31k3p_MASK13) << 7) + (g1[2] >> 24);\n"
   
   "	if (y1 >= mrg31k3p_M1)\n"
   "		y1 -= mrg31k3p_M1;\n"
   
   "	y1 += g1[2];\n"
   "	if (y1 >= mrg31k3p_M1)\n"
   "		y1 -= mrg31k3p_M1;\n"
   
   "	g1[2] = g1[1];\n"
   "	g1[1] = g1[0];\n"
   "	g1[0] = y1;\n"
   
   // second component
   "	y1 = ((g2[0] & mrg31k3p_MASK2) << 15) + (mrg31k3p_MULT2 * (g2[0] >> 16));\n"
   "	if (y1 >= mrg31k3p_M2)\n"
   "		y1 -= mrg31k3p_M2;\n"
   "	y2 = ((g2[2] & mrg31k3p_MASK2) << 15) + (mrg31k3p_MULT2 * (g2[2] >> 16));\n"
   "	if (y2 >= mrg31k3p_M2)\n"
   "		y2 -= mrg31k3p_M2;\n"
   "	y2 += g2[2];\n"
   "	if (y2 >= mrg31k3p_M2)\n"
   "		y2 -= mrg31k3p_M2;\n"
   "	y2 += y1;\n"
   "	if (y2 >= mrg31k3p_M2)\n"
   "		y2 -= mrg31k3p_M2;\n"
   
   "	g2[2] = g2[1];\n"
   "	g2[1] = g2[0];\n"
   "	g2[0] = y2;\n"
   
   "	if (g1[0] <= g2[0]){\n"
   "		temp= g1[0] - g2[0] + mrg31k3p_M1;\n"
   "	} else {\n"
   "		temp = g1[0] - g2[0];\n"
   " }\n";*/
  
  if(random_type == "normal"){  
    result += 
      "    if(get_local_id(1)) {\n"
      "      part[1] = TWOPI_mrg31k3p_NORM_cl * temp;\n"
      "      cosPart1 = cos(part[1]);\n"
      "      sinPart1 = sin(part[1]);\n"
      "    } else {\n"
      "      part[0] = sqrt(-2.0*log(mrg31k3p_NORM_cl * temp));\n"
      "    }\n"
      "    barrier(CLK_LOCAL_MEM_FENCE);\n";
    
    result += 
      "    if(get_local_id(1)) {\n"
      "      out[Dentry] = part[0]*sinPart1;\n"
      "    } else {\n"
      "      out[Dentry] = part[0]*cosPart1;\n"
      "    }\n"
      "    barrier(CLK_LOCAL_MEM_FENCE);\n";
    
  } else if (random_type == "uniform"){ // uniform
    if(typeString == "double" | typeString == "float") {
      result += 
        "out[Dentry] = mrg31k3p_NORM_cl * temp;\n";
    } else if (typeString == "int"){
      result += 
        // "  out[Dentry] = (uint) ((2 * mrg31k3p_M1 + 1) * mrg31k3p_NORM_cl * temp);\n";
        "  out[Dentry] = (int) temp;\n";
    }
    
  }else if(random_type == "exponential"){
    result +=    // "if (mrg31k3p_NORM_cl * temp < 1)\n"
      "out[Dentry] = - log( mrg31k3p_NORM_cl * temp);\n";
  }
  
  result += 
    "  }//Dcol\n";
  result += 
    "}//Drow\n";
  
  result += "streamsFromPrivate(streams,g1,g2,startvalue);\n";
  
  /*
   " for(Drow = 0,DrowStart = index * NpadStreams,Dcol = DrowStart + 3;\n"
   "     Drow < 3; Drow++, DrowStart++, Dcol++){\n"
   "   streams[DrowStart] = g1[Drow];\n"
   "   streams[Dcol] = g2[Drow];\n"
   " }\n";
   */ 
  
  result += 
  "}//kernel\n";
  
  return(result);
}


template<typename T>
int gpuMatrixRn(
    viennacl::matrix<T> &x,
    viennacl::matrix<cl_uint> &streams,
    const IntegerVector numWorkItems,
    const int ctx_id,
    const std::string random_type){
  
  
  std::string mrg31k3pkernelString = mrg31k3pMatrixString<T>(
    x.size1(),
    x.size2(),
    x.internal_size2(),
    streams.internal_size2(),
    random_type);
  
#ifdef DEBUGKERNEL
  Rcpp::Rcout << mrg31k3pkernelString << "\n\n";
#endif  
  
  // the context
  viennacl::ocl::switch_context(ctx_id);
  viennacl::ocl::program & my_prog = viennacl::ocl::current_context().add_program(mrg31k3pkernelString, "my_kernel");
  
  
  viennacl::ocl::kernel &random_number = my_prog.get_kernel("mrg31k3pMatrix");
  
  random_number.global_work_size(0, numWorkItems[0]);
  random_number.global_work_size(1, numWorkItems[1]);
  
  random_number.local_work_size(0, 1L);
  random_number.local_work_size(1, 2L);
  
  
  
  viennacl::ocl::command_queue theQueue = random_number.context().get_queue();
  viennacl::ocl::enqueue(random_number(streams, x), theQueue);
  clFinish(theQueue.handle().get());
  
  return(0L);
}



template<typename T> 
SEXP gpuRnMatrixTyped(
    Rcpp::S4  xR,
    Rcpp::S4  streamsR,
    Rcpp::IntegerVector max_global_size,
    std::string  random_type) 
{
  
  const bool BisVCL=1;
  const int ctx_id = INTEGER(xR.slot(".context_index"))[0]-1;
  
  std::shared_ptr<viennacl::matrix<T> > x = getVCLptr<T>(xR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<cl_uint> > streams = getVCLptr<cl_uint>(streamsR.slot("address"), BisVCL, ctx_id);
  
  
  return(Rcpp::wrap(gpuMatrixRn<T>(*x, *streams, max_global_size, ctx_id, random_type)));	
  
}





// [[Rcpp::export]]
SEXP gpuRnBackend(
    Rcpp::S4  x,
    Rcpp::S4  streams,
    IntegerVector max_global_size,
    std::string  random_type) {
  
  SEXP result;
  
  //Rcpp::Rcout << "222" << "\n\n";
  
  
  Rcpp::traits::input_parameter< std::string >::type classInput(RCPP_GET_CLASS(x));
  std::string classInputString = (std::string) classInput;
  
  
  
  if(classInputString == "fvclMatrix") {
    result = gpuRnMatrixTyped<float>(x, streams, max_global_size, random_type);
  } else if (classInputString == "dvclMatrix") {
    result = gpuRnMatrixTyped<double>(x, streams, max_global_size, random_type);
  } else if (classInputString == "ivclMatrix") {
    result = gpuRnMatrixTyped<int>(x, streams, max_global_size, random_type);
  } else {
    result = Rcpp::wrap(1L);
  }
  return(result);
}






