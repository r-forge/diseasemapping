
std::string mrg31k3pkernelString = 
//"#include <clRNG/clRNG.clh>"

//"#define _CLRNG_FPTYPE cl_double\n"
"#define MODULAR_NUMBER_TYPE cl_uint\n"
"#define CLRNG_ENABLE_SUBSTREAMS\n"


"#define __CLRNG_DEVICE_API\n"
"#pragma OPENCL EXTENSION cl_amd_fp64 : enable\n"

"typedef double cl_double;\n"
"typedef float  cl_float;\n"
"typedef int    cl_int;\n"
"typedef uint   cl_uint;\n"
"typedef long   cl_long;\n"
"typedef ulong  cl_ulong;\n"


//"#define _CLRNG_TAG_FPTYPE(name)           _CLRNG_TAG_FPTYPE_(name,_CLRNG_FPTYPE)\n"
//"#define _CLRNG_TAG_FPTYPE_(name,fptype)   _CLRNG_TAG_FPTYPE__(name,fptype)\n"
//"#define _CLRNG_TAG_FPTYPE__(name,fptype)  name##_##fptype\n"



"typedef enum clrngStatus_ {"
"    CLRNG_SUCCESS              = 0,"
"    CLRNG_INVALID_VALUE        = -1"
"} clrngStatus;\n"

"/*! This macro does nothing."
" *  It is defined for convenience when adapting host code for the device."
" */"
"#define clrngSetErrorString(err) (err)\n"

//"#include <clRNG/mrg31k3p.clh>"

"#define MRG31K3P_CLH\n"

"/********************************************************************************\n"
" * Functions and types declarations                                             *\n"
" ********************************************************************************/\n"

"typedef struct {\n"
"    /*! @brief Seed for the first MRG component\n"
"     */\n"
"    cl_uint g1[3];\n"
"    /*! @brief Seed for the second MRG component\n"
"     */\n"
"    cl_uint g2[3];\n"
"} clrngMrg31k3pStreamState;\n"

"struct clrngMrg31k3pStream_ {\n"
"    clrngMrg31k3pStreamState current;\n"
"#if __OPENCL_C_VERSION__ >= 200\n"
"    // use generic address space\n"
"    const clrngMrg31k3pStreamState* initial;\n"
"#else\n"
"    // force global address space\n"
"    __global const clrngMrg31k3pStreamState* initial;\n"
"#endif\n"
"#ifdef CLRNG_ENABLE_SUBSTREAMS\n"
"    clrngMrg31k3pStreamState substream;\n"
"#endif\n"
"};\n"
"typedef struct clrngMrg31k3pStream_ clrngMrg31k3pStream;\n"


"struct clrngMrg31k3pHostStream_ {\n"
"    clrngMrg31k3pStreamState current;\n"
"    clrngMrg31k3pStreamState initial;\n"
"    clrngMrg31k3pStreamState substream;\n"
"};\n"
"typedef struct clrngMrg31k3pHostStream_ clrngMrg31k3pHostStream;\n"



" /******************************************************************************** \n"
" * Implementation                                                               * \n"
" ********************************************************************************/ \n"



"clrngStatus clrngMrg31k3pCopyOverStreamsFromGlobal(size_t count, clrngMrg31k3pStream* destStreams, __global const clrngMrg31k3pHostStream* srcStreams)\n"
"{\n"
//"    //Check params\n"
//"    if (!destStreams)\n"
//"	return clrngSetErrorString(CLRNG_INVALID_VALUE, \"clrngMrg31k3pCopyOverStreamsFromGlobal(): destStreams cannot be NULL\");\n"
//"    if (!srcStreams)\n"
//"	return clrngSetErrorString(CLRNG_INVALID_VALUE, \"clrngMrg31k3pCopyOverStreamsFromGlobal(): srcStreams cannot be NULL\");\n"

"    for (size_t i = 0; i < count; i++) {\n"
"	destStreams[i].current   = srcStreams[i].current;\n"
"	destStreams[i].initial   = &srcStreams[i].initial;\n"
"#ifdef CLRNG_ENABLE_SUBSTREAMS\n"
"	destStreams[i].substream = srcStreams[i].substream;\n"
"#endif\n"
"    }\n"

"    return CLRNG_SUCCESS;\n"
"}\n"


"clrngStatus clrngMrg31k3pCopyOverStreamsToGlobal(size_t count, __global clrngMrg31k3pHostStream* destStreams, const clrngMrg31k3pStream* srcStreams)\n"
"{\n"
//"    //Check params\n"
//"    if (!destStreams){\n"
//"	return clrngSetErrorString(CLRNG_INVALID_VALUE, \" clrngMrg31k3pCopyOverStreamsToGlobal(): destStreams cannot be NULL \" );\n"
//"}\n"
//"    if (!srcStreams){\n"
//"	return clrngSetErrorString(CLRNG_INVALID_VALUE, \"clrngMrg31k3pCopyOverStreamsToGlobal(): srcStreams cannot be NULL\");\n"
//"}\n"

"    for (size_t i = 0; i < count; i++) {\n"
"	destStreams[i].current   = srcStreams[i].current;\n"
"	destStreams[i].initial   = *srcStreams[i].initial;\n"
"#ifdef CLRNG_ENABLE_SUBSTREAMS\n"
"	destStreams[i].substream = srcStreams[i].substream;\n"
"#endif\n"
"    }\n"

"    return CLRNG_SUCCESS;\n"
"}\n"



// code that is common to host and device
//#include <clRNG/private/mrg31k3p.c.h>

"#define mrg31k3p_M1 2147483647\n"             /* 2^31 - 1 */
"#define mrg31k3p_M2 2147462579\n"             /* 2^31 - 21069 */

"#define mrg31k3p_MASK12 511  \n"              /* 2^9 - 1 */
"#define mrg31k3p_MASK13 16777215  \n"         /* 2^24 - 1 */
"#define mrg31k3p_MASK2 65535     \n"          /* 2^16 - 1 */
"#define mrg31k3p_MULT2 21069\n"

"#define mrg31k3p_NORM_cl_double 4.656612873077392578125e-10 \n" /* 1/2^31 */
"#define mrg31k3p_NORM_cl_float  4.6566126e-10\n"


/*! @brief Advance the rng one step and returns z such that 1 <= z <= mrg31k3p_M1
 */
"static cl_uint clrngMrg31k3pNextState(clrngMrg31k3pStreamState* currentState)\n"
"{\n"
	
"	cl_uint* g1 = currentState->g1;\n"
"	cl_uint* g2 = currentState->g2;\n"
"	cl_uint y1, y2;\n"

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

"	if (g1[0] <= g2[0])\n"
"		return (g1[0] - g2[0] + mrg31k3p_M1);\n"
"	else\n"
"		return (g1[0] - g2[0]);\n"
"}\n"








"__kernel void mrg31k3pDoubleUint("
  "__global clrngMrg31k3pHostStream* streams, __global double* out,\n"
  "const int Nsim) {\n"

"const int Dglobal = get_global_id(0), Dlocal = get_local_id(0);\n"
"const int NlocalSize = get_local_size(0);\n"
"const int Nsize = get_global_size(0);\n"

"clrngMrg31k3pStream private_stream_d;\n" // This is not a pointer!
"clrngMrg31k3pCopyOverStreamsFromGlobal(1, &private_stream_d, &streams[Dglobal]);\n"

"int D;\n"

"for(D = Dglobal; D < Nsim; D += Nsize){\n"
	"out[D] = clrngMrg31k3pNextState(&private_stream_d.current) * mrg31k3p_NORM_cl_double;\n"
"}\n"

"clrngMrg31k3pCopyOverStreamsToGlobal(1,  &streams[Dglobal], &private_stream_d);"
"}\n"






"__kernel void mrg31k3pIntUint("
  "__global clrngMrg31k3pHostStream* streams, __global int* out,\n"
  "const int Nsim) {\n"

"const int Dglobal = get_global_id(0), Dlocal = get_local_id(0);\n"
"const int NlocalSize = get_local_size(0);\n"
"const int Nsize = get_global_size(0);\n"

"clrngMrg31k3pStream private_stream_d;\n" // This is not a pointer!
"clrngMrg31k3pCopyOverStreamsFromGlobal(1, &private_stream_d, &streams[Dglobal]);\n"

"int D;\n"

"for(D = Dglobal; D < Nsim; D += Nsize){\n"
	"out[D] = (int) clrngMrg31k3pNextState(&private_stream_d.current);\n"
"}\n"

"clrngMrg31k3pCopyOverStreamsToGlobal(1,  &streams[Dglobal], &private_stream_d);"
"}"









"__kernel void mrg31k3pFloatUint("
  "__global clrngMrg31k3pHostStream* streams, __global float* out,\n"
  "const int Nsim) {\n"

"const int Dglobal = get_global_id(0), Dlocal = get_local_id(0);\n"
"const int NlocalSize = get_local_size(0);\n"
"const int Nsize = get_global_size(0);\n"

"clrngMrg31k3pStream private_stream_d;\n" // This is not a pointer! the declaration allocates private memory
"clrngMrg31k3pCopyOverStreamsFromGlobal(1, &private_stream_d, &streams[Dglobal]);\n" //copy from host into private memory

"int D;\n"

"for(D = Dglobal; D < Nsim; D += Nsize){\n"
	"out[D] = clrngMrg31k3pNextState(&private_stream_d.current) * mrg31k3p_NORM_cl_float;\n"
"}\n"

"clrngMrg31k3pCopyOverStreamsToGlobal(1,  &streams[Dglobal], &private_stream_d);"//copy from device private into global 
"}";




  
  "__kernel void mrg31k3pFloatNorm("
  "__global clrngMrg31k3pHostStream* streams, __global float* u1,\n"
  "const int Nsim) {\n"
  
  "const int Dglobal = get_global_id(0);\n"
  "const int Nsize = get_global_size(0);\n"
  
  "clrngMrg31k3pStream private_stream_d;\n" // This is not a pointer! the declaration allocates private memory
  "clrngMrg31k3pCopyOverStreamsFromGlobal(1, &private_stream_d, &streams[Dglobal]);\n" //copy from host into private memory
  
  "int i;\n"
  
  "for(i = Dglobal; i < Nsim; i += Nsize){\n"
  "u1[i] = clrngMrg31k3pRandomU01(&private_stream_d);\n"
  "}\n"
  
  "clrngMrg31k3pCopyOverStreamsToGlobal(1,  &streams[Dglobal], &private_stream_d);"//copy from device private into global 
  "}";





