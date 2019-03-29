
std::string mrg31k3pPrefixDouble = 
"#pragma OPENCL EXTENSION cl_amd_fp64 : enable\n"
"#define mrg31k3p_convert(xx) ( (xx) * 4.656612873077392578125e-10)\n" /* 1/2^31 */
"#define OUTPUT_TYPE double\n"
"\n";


std::string mrg31k3pPrefixFloat = 
"#define mrg31k3p_convert(xx) ( (xx) * 4.6566126e-10)\n"
"#define OUTPUT_TYPE float\n"
"\n";

std::string mrg31k3pPrefixInt = 
"#define OUTPUT_TYPE int\n"
"#define mrg31k3p_convert(xx) ( (int) xx )\n"
"\n";



std::string mrg31k3pTemplateString = 

"\n\n__kernel void mrg31k3p("
  "__global clrngMrg31k3pHostStream* streams," 
  "__global OUTPUT_TYPE* out,\n"
  "const int Nsim) {\n"

"const int Dglobal = get_global_id(0);\n" 
"const int Nsize = get_global_size(0);\n"

"clrngMrg31k3pStream private_stream_d;\n" // This is not a pointer!
"clrngMrg31k3pCopyOverStreamsFromGlobal(1, &private_stream_d, &streams[Dglobal]);\n"

"int D;\n"

"for(D = Dglobal; D < Nsim; D += Nsize){\n"
	"out[D] = mrg31k3p_convert(clrngMrg31k3pNextState(&private_stream_d.current));\n"
"}\n"

"clrngMrg31k3pCopyOverStreamsToGlobal(1,  &streams[Dglobal], &private_stream_d);"
"}\n"
"\n";

