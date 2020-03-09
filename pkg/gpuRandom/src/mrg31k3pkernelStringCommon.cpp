#include <string>

std::string mrg31k3pString() {

std::string result = "\n#define CLRNG_ENABLE_SUBSTREAMS\n"
"#define __CLRNG_DEVICE_API\n"

"#define mrg31k3p_M1 2147483647\n"             /* 2^31 - 1 */
"#define mrg31k3p_M2 2147462579\n"             /* 2^31 - 21069 */

"#define mrg31k3p_MASK12 511  \n"              /* 2^9 - 1 */
"#define mrg31k3p_MASK13 16777215  \n"         /* 2^24 - 1 */
"#define mrg31k3p_MASK2 65535     \n"          /* 2^16 - 1 */
"#define mrg31k3p_MULT2 21069\n"

//"#define MODULAR_NUMBER_TYPE cl_uint\n"

//"typedef double cl_double;\n"
//"typedef float  cl_float;\n"
//"typedef int    cl_int;\n"
// "typedef uint   cl_uint;\n"
//"typedef long   cl_long;\n"
//"typedef ulong  cl_ulong;\n"


//"#define MRG31K3P_CLH\n"


" /******************************************************************************** \n"
" * Implementation                                                               * \n"
" ********************************************************************************/ \n"

"void streamsToPrivate(__global int* streams, uint* g1, uint* g2, const int start){\n"
" int Drow, Dcol, DrowStart;"
"for(Drow = 0, DrowStart = start, Dcol = DrowStart + 3;\n"
"Drow < 3; Drow++, DrowStart++, Dcol++){\n"
"g1[Drow] = streams[DrowStart];\n"
"g2[Drow] = streams[Dcol];\n"
"}\n"
"}\n"


"void streamsFromPrivate(__global int* streams, uint* g1, uint* g2,  const int start){\n"
" int Drow, Dcol, DrowStart;"
" for(Drow = 0,DrowStart = start, Dcol = DrowStart + 3;\n"
"     Drow < 3; Drow++, DrowStart++, Dcol++){\n"
"   streams[DrowStart] = g1[Drow];\n"
"   streams[Dcol] = g2[Drow];\n"
" }\n"
"}\n"



/*! @brief Advance the rng one step and returns z such that 1 <= z <= mrg31k3p_M1
 */


"uint clrngMrg31k3pNextState(uint *g1, uint *g2) {\n"

"	uint y1, y2;\n"

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
"		return (g1[0] - g2[0] + mrg31k3p_M1);\n"      //"		result = convert_int(mrg31k3p_M1 - (g2[0] - g1[0]));\n"
"	} else {\n"
"		return(g1[0] - g2[0]);\n"
" }\n"
"}\n\n";

  

  
  
  
return(result);
}







