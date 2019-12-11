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
    "\n#define mrg31k3p_NORM_cl 4.6566126e-10\n\n";
    "//TWOPI * mrg31k3p_NORM_cl\n"
    "#define TWOPI_mrg31k3p_NORM_cl 2.9258361e-09\n\n";
  } else {
    result += 
    "\n#define mrg31k3p_NORM_cl 1L\n\n";
  }
  
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

  result += 
    "int Drow, Dcol, DrowStart, Dentry;\n";
  
  result += "cl_uint g1[3], g2[3];\n"
  " for(Drow = 0,DrowStart = index * NpadStreams,Dcol = DrowStart + 3;\n"
  "     Drow < 3; Drow++, DrowStart++, Dcol++){"
  "   g1[Drow] = streams[DrowStart];\n"
  "   g2[Drow] = streams[Dcol];\n"
  " }\n";    
    result += "cl_uint temp;\n";
  
    if(random_type == "normal"){  
      result += 
      "local " + typeString + " part[2];\n";// local size must be 1,2
    }
  
  result += "const int DrowInc = 2*get_global_size(0),\n"
  "  DrowStartInc = DrowInc * NpadCol;\n";
  
  
  // Drow, Dcol = 2* global(0), group(1)
  // local work item (0,0) does out[Drow, Dcol]
  // local work item (0,1) does out[Drow+1, Dcol]
  result += 
    "for(Drow=2*get_global_id(0), DrowStart = (Drow + get_local_id(1))* NpadCol;\n" 
    "    Drow < Nrow; Drow += DrowInc, DrowStart += DrowStartInc) {\n";
  
  result += 
    "  for(Dcol=get_group_id(1), Dentry = DrowStart + Dcol;\n" 
    "      Dcol < Ncol;\n" 
    "      Dcol += get_num_groups(1), Dentry += get_num_groups(1)) {\n";
  
  result += 
    "    temp = clrngMrg31k3pNextState(g1, g2);\n";

  if(random_type == "normal"){  
  result += 
    "    if(get_local_id(1)) {\n"
    "      part[1] = TWOPI_mrg31k3p_NORM_cl * temp;\n"
    //"sinPart1 = sin(part[1]);\n"
    //"part[1] = cos(part[1]);\n"
    "    } else {"
    "      part[0] = sqrt(-2.0*log(temp * mrg31k3p_NORM_cl));\n"
    "    }\n"
    "    barrier(CLK_LOCAL_MEM_FENCE);\n";
  
  result += 
    "    if(get_local_id(1)) {\n"
    "      out[Dentry] = part[0]*sin(part[1]);\n"
    "    } else {"
    "      out[Dentry] = part[0]*cos(part[1]);\n"
    "    }\n"
    "    barrier(CLK_LOCAL_MEM_FENCE);\n";
  } else { // uniform
    if(typeString == "double" | typeString == "float") {
    result += 
      "out[Dentry] = mrg31k3p_NORM_cl * temp;\n";
    } else {
      result += 
        "  out[Dentry] = temp;\n";
    }
  }
  
  result += 
    "  }//Dcol\n "
    "}//Drow\n";
  
  result +=
    " for(Drow = 0,DrowStart = index * NpadStreams,Dcol = DrowStart + 3;\n"
    "     Drow < 3; Drow++, DrowStart++, Dcol++){"
    "   streams[DrowStart] = g1[Drow];\n"
    "   streams[Dcol] = g2[Drow];\n"
    " }\n";
  
  result += 
    "}//kernel\n";
  
  return(result);
}


template<typename T>
void gpuMatrixRn(
    viennacl::matrix<T> &x,
    const Rcpp::IntegerMatrix streamsIn,
    Rcpp::IntegerMatrix streamsOut,
    const IntegerVector numWorkItems,
    const int ctx_id,
    const std::string random_type){

   const int Nstreams = numWorkItems[0]*numWorkItems[1];
   viennacl::matrix<int> streamsGpu(Nstreams, 6);
   int Drow, Dcol;

   for(Drow = 0; Drow < Nstreams; Drow++) {
     for(Dcol = 0; Dcol < 6; Dcol++) {
        streamsGpu(Drow, Dcol) = streamsIn(Drow, Dcol);       
     }
   }

  std::string mrg31k3pkernelString = mrg31k3pMatrixString<T>(
    x.size1(),
    x.size2(),
    x.internal_size2(),
    streamsGpu.internal_size2(),
    random_type);
  
#ifdef DEBUGKERNEL
  
  Rcpp::Rcout << mrg31k3pkernelString << "\n\n";
  
#endif  

  
  viennacl::ocl::context ctx(viennacl::ocl::get_context(ctx_id));
  viennacl::ocl::program &my_prog = ctx.add_program(mrg31k3pkernelString, "my_kernel");
  
  viennacl::ocl::kernel &random_number = my_prog.get_kernel("mrg31k3pMatrix");
  random_number.global_work_size(0, numWorkItems[0]);
  random_number.global_work_size(1, numWorkItems[1]);
  
  random_number.local_work_size(0, 1L);
  random_number.local_work_size(1, 2L);

  viennacl::ocl::enqueue(random_number(streamsGpu, x) );


#ifdef DEBUG
  
  Rcpp::Rcout << "x\n" << x(0,0) << " "<< x(0,1) << " // "<< x(1,0) << " "<< x(1,1) <<  "\n";
//  Rcpp::Rcout << x2.size1() << " "<< x2.size2() << "  "<< x2.internal_size1() << " "<< x.internal_size2()  <<  "\n";
  
#endif  
  

  // copy streams back to cpu
  for(Drow = 0; Drow < Nstreams; Drow++) {
    for(Dcol = 0; Dcol < 6; Dcol++) {
      streamsOut(Drow, Dcol) = streamsGpu(Drow, Dcol);       
    }
    for(Dcol = 6 ; Dcol < 18; Dcol++) {
      streamsOut(Drow, Dcol) = streamsIn(Drow, Dcol);       
    }
  }

}



template<typename T> 
SEXP gpuRnMatrixTyped(
    Rcpp::S4  xR,
    const Rcpp::IntegerMatrix streamsIn,
    IntegerVector max_global_size,
    std::string  random_type) 
{
  
  Rcpp::IntegerMatrix streamsOut(streamsIn.nrow(), streamsIn.ncol());

  const bool BisVCL=1;
  const int ctx_id = INTEGER(xR.slot(".context_index"))[0]-1;
  
  std::shared_ptr<viennacl::matrix<T> > x = getVCLptr<T>(
    xR.slot("address"), BisVCL, ctx_id);
  
  gpuMatrixRn<T>(*x, streamsIn, streamsOut, 
                 max_global_size, ctx_id, random_type);
  
  return(streamsOut);	
}


//' Random number generation
//' 
//' Fills a matrix with random numbers
//'
//' @param x output matrix
//' @param streams matrix of random number seeds
//' @param max_global_size vector of length 2, number of work items
//' @param random_type one of "uniform" or "normal"
//' 
//' @export
// [[Rcpp::export]]
SEXP gpuRnBackend(
    Rcpp::S4  x,
    const Rcpp::IntegerMatrix streams,
    IntegerVector max_global_size,
    std::string  random_type) {
  
  SEXP result;
  
  Rcpp::traits::input_parameter< std::string >::type classInput(RCPP_GET_CLASS(x));
  std::string classInputString = (std::string) classInput;
  
  
  if(classInputString == "fvclMatrix") {
    result = gpuRnMatrixTyped<float>(
      x, streams, max_global_size, random_type);
  } else if (classInputString == "dvclMatrix") {
    result = gpuRnMatrixTyped<double>(
      x, streams, max_global_size, random_type);
  } else if (classInputString == "ivclMatrix") {
    result = gpuRnMatrixTyped<int>(
      x, streams, max_global_size, random_type);
  } else {
    result = Rcpp::wrap(1L);
  }
  return(result);
}
