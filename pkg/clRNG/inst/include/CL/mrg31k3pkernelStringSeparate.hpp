


std::string mrg31k3pTemplateStringFirst = 

"\n\n__kernel void mrg31k3p(\n"
"  __global clrngMrg31k3pHostStream* streams,\n" 
"  __global "; //OUTPUT_TYPE



std::string mrg31k3pTemplateStringSecond = 
"* out,\n"
"const int Nsim) {\n"

"const int Xsize = get_global_size(0);\n"
"const int size = get_global_size(1)*Xsize;\n"

"clrngMrg31k3pStream private_stream_d;\n" // This is not a pointer! the declaration allocates private memory
"clrngMrg31k3pCopyOverStreamsFromGlobal(1, &private_stream_d, &streams[Dglobal]);\n" //copy from host into private memory

"int index = get_global_id(1)*Xsize + get_global_id(0);\n"
"int D;\n"

"for (D=index; D < Nsim ; D += size) {\n"
"out[D] ="



std::string mrg31k3pTemplateStringThird = 
"clrngMrg31k3pNextState(&private_stream_d.current) ;\n"
"}\n"

"clrngMrg31k3pCopyOverStreamsToGlobal(1,  &streams[Dglobal], &private_stream_d);"//Copy RNG device stream objects from private memory into global memory
"}\n"


  // normal float string and normal double string
std::string mrg31k3pTemplateStringForth = 
"* z,\n"
"const int Nsim ){\n"

"const int Xsize = get_global_size(0);\n"
"const int size = get_global_size(1)*Xsize;\n"

"clrngMrg31k3pStream private_stream_d;\n"// This is not a pointer! the declaration allocates private memory
"clrngMrg31k3pCopyOverStreamsFromGlobal(1, &private_stream_d, &streams[Dglobal]);\n"//copy from host into private memory

"int index = (get_global_id(1)*Xsize + get_global_id(0))*2\n"
"int i, D;\n"

// multiply by the relevant constant

std::string mrg31k3pTemplateStringFifth = 
"temp[2];\n"
"temp[0] = clrngMrg31k3pNextState(&private_stream_d.current) * mrg31k3p_NORM_cl;\n "//clrngMrg31k3pRandomU01(&private_stream_d);\n"
"temp[1] = clrngMrg31k3pNextState(&private_stream_d.current) * mrg31k3p_NORM_cl;\n"

"for (D=index; D < Nsim; D+=2*size) {\n"
   "for (i = 0; i < 2; ++i) {\n"
      "if (D+i < Nsim) {\n"
      "z[D+i] = sqrt(-2.0f*log(temp[0]))*((1-i)*cos(TWOPI*temp[1])+i*sin(TWOPI*temp[1]));\n"
 "}\n"
"}\n"
"}\n"
 
"clrngMrg31k3pCopyOverStreamsToGlobal(1,  &streams[Dglobal], &private_stream_d);\n" //a single stream object
"};\n"
;




std::string mrg31k3pDoubleUnifString = 
    "\n#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n\n" + 
    mrg31k3pCommon +
    mrg31k3pTemplateStringFirst + "double" +
    mrg31k3pTemplateStringSecond +
    "4.656612873077392578125e-10 * " +
    mrg31k3pTemplateStringThird;


std::string mrg31k3pFloatUnifString = 
    "\n#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n\n" + 
    mrg31k3pCommon + 
    mrg31k3pTemplateStringFirst + "float" +
    mrg31k3pTemplateStringSecond +
    "4.6566126e-10 * " +
    mrg31k3pTemplateStringThird;


std::string mrg31k3pIntegerUnifString = 
    mrg31k3pCommon + 
    mrg31k3pTemplateStringFirst + "int" +
    mrg31k3pTemplateStringSecond + 
    "(int) " + 
    mrg31k3pTemplateStringThird;



///////////////////////////////////////////////////////////////////////
std::string mrg31k3pDoubleNormString = 
    "\n#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n\n" + 
    "\n#define TWOPI 6.283185307179586 \n\n" + 
    "\n#define mrg31k3p_NORM_cl 4.656612873077392578125e-10\n\n"+
    mrg31k3pCommon +
    mrg31k3pTemplateStringFirst + "double" +
    mrg31k3pTemplateStringForth + "double" +
    mrg31k3pTemplateStringFifth;


std::string mrg31k3pFloatNormString = 
    "\n#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n\n" + 
    "\n#define TWOPI 6.2831853\n\n" + 
    "\n#define mrg31k3p_NORM_cl 4.6566126e-10\n\n"+
    mrg31k3pCommon +
    mrg31k3pTemplateStringFirst + "float" +
    mrg31k3pTemplateStringForth + "float" +
    mrg31k3pTemplateStringFifth;









