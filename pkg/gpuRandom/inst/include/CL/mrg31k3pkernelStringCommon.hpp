#ifndef mrg31k3pCommonDefined

#define mrg31k3pCommonDefined

std::string mrg31k3pCommon =
"#define CLRNG_ENABLE_SUBSTREAMS\n"
"#define __CLRNG_DEVICE_API\n"

"#define mrg31k3p_M1 2147483647\n"             /* 2^31 - 1 */
"#define mrg31k3p_M2 2147462579\n"             /* 2^31 - 21069 */

"#define mrg31k3p_MASK12 511  \n"              /* 2^9 - 1 */
"#define mrg31k3p_MASK13 16777215  \n"         /* 2^24 - 1 */
"#define mrg31k3p_MASK2 65535     \n"          /* 2^16 - 1 */
"#define mrg31k3p_MULT2 21069\n"

"#define MODULAR_NUMBER_TYPE cl_uint\n"

//"typedef double cl_double;\n"
"typedef float  cl_float;\n"
"typedef int    cl_int;\n"
"typedef uint   cl_uint;\n"
"typedef long   cl_long;\n"
"typedef ulong  cl_ulong;\n"


"typedef enum clrngStatus_ {\n"
"    CLRNG_SUCCESS              = 0,\n"
"    CLRNG_INVALID_VALUE        = -1\n"
"} clrngStatus;\n"

"/*! This macro does nothing."
" *  It is defined for convenience when adapting host code for the device."
" */\n"
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

"    for (size_t i = 0; i < count; i++) {\n"
"	destStreams[i].current   = srcStreams[i].current;\n"
"	destStreams[i].initial   = *srcStreams[i].initial;\n"
"#ifdef CLRNG_ENABLE_SUBSTREAMS\n"
"	destStreams[i].substream = srcStreams[i].substream;\n"
"#endif\n"
"    }\n"

"    return CLRNG_SUCCESS;\n"
"}\n"
/*! @brief Advance the rng one step and returns z such that 1 <= z <= mrg31k3p_M1
 */
//"static 
"cl_uint clrngMrg31k3pNextState(clrngMrg31k3pStreamState* currentState) {\n"
	
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
"\n";

#endif