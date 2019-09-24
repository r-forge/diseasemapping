#include <CL/mrg31k3pkernelStringCommon.hpp>


std::string mrg31k3pTemplateStringFirst = 

"\n\n__kernel void mrg31k3p(\n"
"  __global clrngMrg31k3pHostStream* streams,\n" 
"  __global "; //OUTPUT_TYPE



std::string mrg31k3pTemplateStringSecond = 
"* out,\n"
"const int Nsim) {\n"

"const int size = get_global_size(1)*get_global_size(0);\n"
"int index = get_global_id(1)*get_global_size(0) + get_global_id(0);\n"

"clrngMrg31k3pStream private_stream_d;\n" // This is not a pointer! the declaration allocates private memory
"clrngMrg31k3pCopyOverStreamsFromGlobal(1, &private_stream_d, &streams[index]);\n" //copy from host into private memory

"int D;\n"

"for (D=index; D < Nsim ; D += size) {\n"
"out[D] = ";



std::string mrg31k3pTemplateStringThird = 
"clrngMrg31k3pNextState(&private_stream_d.current) ;\n"
"}\n"

"clrngMrg31k3pCopyOverStreamsToGlobal(1,  &streams[index], &private_stream_d);"//Copy RNG device stream objects from private memory into global memory
"}\n";


  // normal float string and normal double string
std::string mrg31k3pTemplateStringForth = 
"* z,\n"
"const int Nsim ){\n"

"const int size = get_global_size(1)*get_global_size(0);\n"

"clrngMrg31k3pStream private_stream_d;\n"// This is not a pointer! the declaration allocates private memory
"int index = get_global_id(1)*get_global_size(0) + get_global_id(0);\n"
"int D;\n"
"clrngMrg31k3pCopyOverStreamsFromGlobal(1, &private_stream_d, &streams[index]);\n"
"index = 2*index;\n";
//copy from host into private memory



std::string mrg31k3pTemplateStringFifth = 
" "
"temp0,temp1,part1,part2;\n"
"for (D=index; D < Nsim; D+=2*size) {\n"

  "temp0 = clrngMrg31k3pNextState(&private_stream_d.current)* mrg31k3p_NORM_cl;\n"    //"clrngMrg31k3pRandomU01(&private_stream_d);\n"
  "temp1 = clrngMrg31k3pNextState(&private_stream_d.current)* mrg31k3p_NORM_cl;\n"     //"clrngMrg31k3pRandomU01(&private_stream_d);\n"
  
  "part1 = sqrt(-2.0*log(temp0));\n"
  "part2 = TWOPI * temp1;\n"
  
  "if (D<Nsim) {\n"
   "z[D] = part1*cos(part2);\n"
  "}\n"
   
  "if (D +1< Nsim) {\n"
     "z[D+1] = part1*sin(part2);\n"
  "}\n"
 
 "}\n"
 
"clrngMrg31k3pCopyOverStreamsToGlobal(1,  &streams[index], &private_stream_d);\n" //a single stream object
"}\n"
;




std::string mrg31k3pDoubleUnifString = 
    "\n#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n"
    "\n#define mrg31k3p_NORM_cl 4.656612873077392578125e-10\n\n" +
    mrg31k3pCommon +
    mrg31k3pTemplateStringFirst + "double" +
    mrg31k3pTemplateStringSecond +
    "mrg31k3p_NORM_cl * " +
    mrg31k3pTemplateStringThird;


std::string mrg31k3pFloatUnifString = 
  "\n#define mrg31k3p_NORM_cl 4.6566126e-10\n\n" +
    mrg31k3pCommon + 
    mrg31k3pTemplateStringFirst + "float" +
    mrg31k3pTemplateStringSecond +
    "mrg31k3p_NORM_cl * " +
    mrg31k3pTemplateStringThird;


std::string mrg31k3pIntegerUnifString = 
    mrg31k3pCommon + 
    mrg31k3pTemplateStringFirst + "int" +
    mrg31k3pTemplateStringSecond + 
    "(int) " + 
    mrg31k3pTemplateStringThird;


//////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
std::string mrg31k3pDoubleNormString = 
    "\n#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n\n" 
    "\n#define TWOPI 6.283185307179586 \n\n" 
    "\n#define mrg31k3p_NORM_cl 4.656612873077392578125e-10\n\n"+
    mrg31k3pCommon +
    mrg31k3pTemplateStringFirst + "double" +
    mrg31k3pTemplateStringForth + "double" +
    mrg31k3pTemplateStringFifth;


std::string mrg31k3pFloatNormString = 
    "\n#define TWOPI 6.2831853\n\n" 
    "\n#define mrg31k3p_NORM_cl 4.6566126e-10\n\n"+
    mrg31k3pCommon +
    mrg31k3pTemplateStringFirst + "float" +
    mrg31k3pTemplateStringForth + "float" +
    mrg31k3pTemplateStringFifth;









