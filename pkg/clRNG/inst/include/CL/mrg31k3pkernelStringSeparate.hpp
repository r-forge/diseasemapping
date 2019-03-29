


std::string mrg31k3pTemplateStringFirst = 

"\n\n__kernel void mrg31k3p(\n"
"  __global clrngMrg31k3pHostStream* streams,\n" 
"  __global "; //OUTPUT_TYPE

std::string mrg31k3pTemplateStringSecond = 
"* out,\n"
"  const int Nsim) {\n"

"const int Dglobal = get_global_id(0);\n" 
"const int Nsize = get_global_size(0);\n"

"clrngMrg31k3pStream private_stream_d;\n" // This is not a pointer!
"clrngMrg31k3pCopyOverStreamsFromGlobal(1, &private_stream_d, &streams[Dglobal]);\n"

"int D;\n"

"for(D = Dglobal; D < Nsim; D += Nsize){\n"
	"out[D] = ";

// multiply by the relevant constant

std::string mrg31k3pTemplateStringThird = 
  "clrngMrg31k3pNextState(&private_stream_d.current);\n}\n"
"clrngMrg31k3pCopyOverStreamsToGlobal(1,  &streams[Dglobal], &private_stream_d);\n"
"}\n"
"\n";

