#ifdef UNDEF

//"#include "./private/modular.c.h"\n"

//"#ifndef MODULAR_NUMBER_TYPE\n"
//"#error \"MODULAR_NUMBER_TYPE must be defined\"\n"
//"#endif\n"

//"#ifdef MODULAR_FIXED_SIZE\n"
"#define N MODULAR_FIXED_SIZE\n"
//"#define MATRIX_ELEM(mat, i, j) (mat[i][j])\n"
//"#else\n"
"#define MATRIX_ELEM(mat, i, j) (mat[i * MODULAR_FIXED_SIZE + j])\n"
//"#endif // MODULAR_FIXED_SIZE\n"
"//! Compute (a*s + c) % m\n"
"#if 1\n"
"#define modMult(a, s, c, m) ((MODULAR_NUMBER_TYPE)(((cl_ulong) a * s + c) % m))\n"
"#else\n"
"static MODULAR_NUMBER_TYPE modMult(MODULAR_NUMBER_TYPE a, MODULAR_NUMBER_TYPE s, MODULAR_NUMBER_TYPE c, MODULAR_NUMBER_TYPE m)\n"
"{\n"
"    MODULAR_NUMBER_TYPE v;\n"
"    v = (MODULAR_NUMBER_TYPE) (((cl_ulong) a * s + c) % m);\n"
"    return v;\n"
"}\n"
"#endif\n"



"//! @brief Matrix-vector modular multiplication\n"
"//  @details Also works if v = s.\n"
"//  @return v = A*s % m\n"
//"#ifdef MODULAR_FIXED_SIZE\n"
//"static 
//"void modMatVec(MODULAR_NUMBER_TYPE A[N][N]), MODULAR_NUMBER_TYPE s[N], MODULAR_NUMBER_TYPE v[N], MODULAR_NUMBER_TYPE m)\n"
//"#else\n"
//"void modMatVec(size_t N, MODULAR_NUMBER_TYPE* A)"//, MODULAR_NUMBER_TYPE* s, MODULAR_NUMBER_TYPE* v, MODULAR_NUMBER_TYPE m)\n"
//"#endif\n"
"void modMatVec(local cl_uint *A,"
"	MODULAR_NUMBER_TYPE* s, MODULAR_NUMBER_TYPE* v, MODULAR_NUMBER_TYPE m) {\n"
"   MODULAR_NUMBER_TYPE x[MODULAR_FIXED_SIZE];\n"     // Necessary if v = s
"   for (size_t i = 0; i < MODULAR_FIXED_SIZE; ++i) {\n"
"        x[i] = 0;\n"
"        for (size_t j = 0; j < MODULAR_FIXED_SIZE; j++)\n"
"            x[i] = modMult(MATRIX_ELEM(A,i,j), s[j], x[i], m);\n"
"   }\n"
"    for (size_t i = 0; i < MODULAR_FIXED_SIZE; ++i)\n"
"        v[i] = x[i];\n"
"}\n"

"#undef MATRIX_ELEM\n"
"#undef N\n"

#endif



#ifdef UNDEF
"clrngStatus clrngMrg31k3pCopyOverStreams(size_t count, clrngMrg31k3pStream* destStreams, const clrngMrg31k3pStream* srcStreams)\n"
"{\n"
    //Check params\n"
//	if (!destStreams)\n"
//	return clrngSetErrorString(CLRNG_INVALID_VALUE, "%s(): destStreams cannot be NULL", __func__);\n"
//	if (!srcStreams)\n"
//	return clrngSetErrorString(CLRNG_INVALID_VALUE, "%s(): srcStreams cannot be NULL", __func__);\n"

 "   for (size_t i = 0; i < count; i++)\n"
"		destStreams[i] = srcStreams[i];\n"

"    return CLRNG_SUCCESS;\n"
"}\n"
#endif

std::vector<cl_uint> clrngMrg31k3p_All = { 
    1516919229,  758510237, 499121365,
    1884998244, 1516919229, 335398200,
    601897748,  1884998244, 358115744,
    1228857673, 1496414766,  954677935,
    1133297478, 1407477216, 1496414766,
    2002613992, 1639496704, 1407477216};




#ifdef UNDEF

"clrngStatus clrngMrg31k3pRewindStreams(size_t count, clrngMrg31k3pStream* streams) \n"
"{ \n"
	//Check params
//"	if (!streams) \n"
//"		return clrngSetErrorString(CLRNG_INVALID_VALUE, "%s(): streams cannot be NULL", __func__); \n"
	//Reset current state to the stream initial state
"	for (size_t j = 0; j < count; j++) { \n"
"#ifdef __CLRNG_DEVICE_API \n"
"#ifdef CLRNG_ENABLE_SUBSTREAMS \n"
"		streams[j].current = streams[j].substream = *streams[j].initial; \n"
"#else \n"
"		streams[j].current = *streams[j].initial; \n"
"#endif \n"
"#else \n"
"		streams[j].current = streams[j].substream = streams[j].initial; \n"
"#endif \n"
"	} \n"

"	return CLRNG_SUCCESS; \n"
"} \n"

"#if defined(CLRNG_ENABLE_SUBSTREAMS) || !defined(__CLRNG_DEVICE_API) \n"
"clrngStatus clrngMrg31k3pRewindSubstreams(size_t count, clrngMrg31k3pStream* streams) \n"
"{ \n"
	//Check params
//"	if (!streams) \n"
//"		return clrngSetErrorString(CLRNG_INVALID_VALUE, "%s(): streams cannot be NULL", __func__); \n"
	//Reset current state to the subStream initial state \n"
"	for (size_t j = 0; j < count; j++) { \n"
"		streams[j].current = streams[j].substream; \n"
"	} \n"

"	return CLRNG_SUCCESS; \n"
"} \n"



"clrngStatus clrngMrg31k3pForwardToNextSubstreams(size_t count, clrngMrg31k3pStream* streams,"
"__local MODULAR_NUMBER_TYPE *seedConstants) \n"
"{ \n"
	//Check params
//	if (!streams)
//		return clrngSetErrorString(CLRNG_INVALID_VALUE, "%s(): streams cannot be NULL", __func__);

	
"	for (size_t k = 0; k < count; k++) {"
//"		modMatVec(seedConstants, streams[k].substream.g1, streams[k].substream.g1, mrg31k3p_M1); \n"
//"		modMatVec(clrngMrg31k3p_A2p72, streams[k].substream.g2, streams[k].substream.g2, mrg31k3p_M2); \n"
"		streams[k].current = streams[k].substream; \n"

"	} \n"

"	return CLRNG_SUCCESS; \n"
"} \n"


"clrngStatus clrngMrg31k3pMakeOverSubstreams(clrngMrg31k3pStream* stream, size_t count, clrngMrg31k3pStream* substreams) \n"
"{ \n"
"	for (size_t i = 0; i < count; i++) { \n"
"		clrngStatus err; \n"
		// snapshot current stream into substreams[i]
"		err = clrngMrg31k3pCopyOverStreams(1, &substreams[i], stream); \n"
"		if (err != CLRNG_SUCCESS) \n"
"		    return err; \n"
		// advance to next substream
"		err = clrngMrg31k3pForwardToNextSubstreams(1, stream); \n"
"		if (err != CLRNG_SUCCESS) \n"
"		    return err; \n"
"	}\n"
"	return CLRNG_SUCCESS; \n"
"} \n"

"#endif  \n"// substreams
#endif


// The following would be much cleaner with C++ templates instead of macros.

// We use an underscore on the r.h.s. to avoid potential recursion with certain
// preprocessors.

#ifdef USE_ORIGINAL_MACROS_SHOULD_BE_FALSE
"#define IMPLEMENT_GENERATE_FOR_TYPE(fptype) \\ \n"
"	\\ \n"
"	fptype clrngMrg31k3pRandomU01_##fptype(clrngMrg31k3pStream* stream, __local MODULAR_NUMBER_TYPE *seedConstants) { \\ \n"
//"	    return clrngMrg31k3pNextState(&stream->current, seedConstants) * mrg31k3p_NORM_##fptype; \\ \n"
"	} \\ \n"
"	\\ \n"
"	cl_int clrngMrg31k3pRandomInteger_##fptype(clrngMrg31k3pStream* stream, cl_int i, cl_int j) { \\ \n"
"	    return i + (cl_int)((j - i + 1) * clrngMrg31k3pRandomU01_##fptype(stream)); \\ \n"
"	} \\ \n"
"	\\ \n"
"	clrngStatus clrngMrg31k3pRandomU01Array_##fptype(clrngMrg31k3pStream* stream, size_t count, fptype* buffer) { \\ \n"
//"		if (!stream) \ \n"
//"			return clrngSetErrorString(CLRNG_INVALID_VALUE, "%s(): stream cannot be NULL", __func__); \ \n"
//"		if (!buffer) \ \n"
//"			return clrngSetErrorString(CLRNG_INVALID_VALUE, "%s(): buffer cannot be NULL", __func__); \ \n"
"		for (size_t i = 0; i < count; i++)  \\ \n"
"			buffer[i] = clrngMrg31k3pRandomU01_##fptype(stream); \\ \n"
"		return CLRNG_SUCCESS; \\ \n"
"	} \\ \n"
"	\\ \n"
"	clrngStatus clrngMrg31k3pRandomIntegerArray_##fptype(clrngMrg31k3pStream* stream, cl_int i, cl_int j, size_t count, cl_int* buffer) { \\ \n"
//"		if (!stream) \ \n"
//"			return clrngSetErrorString(CLRNG_INVALID_VALUE, "%s(): stream cannot be NULL", __func__); \ \n"
//"		if (!buffer) \ \n"
//"			return clrngSetErrorString(CLRNG_INVALID_VALUE, "%s(): buffer cannot be NULL", __func__); \ \n"
"		for (size_t k = 0; k < count; k++) \\ \n"
"			buffer[k] = clrngMrg31k3pRandomInteger_##fptype(stream, i, j); \\ \n"
"		return CLRNG_SUCCESS; \\ \n"
"	} \n"

// On the host, implement everything.
// On the device, implement only what is required to avoid cluttering memory.
"#if defined(CLRNG_SINGLE_PRECISION)  || !defined(__CLRNG_DEVICE_API) \n"
"IMPLEMENT_GENERATE_FOR_TYPE(cl_float) \n"
"#endif \n"
"#if !defined(CLRNG_SINGLE_PRECISION) || !defined(__CLRNG_DEVICE_API) \n"
"IMPLEMENT_GENERATE_FOR_TYPE(cl_double) \n"
"#endif \n"

// Clean up macros, especially to avoid polluting device code.
"#undef IMPLEMENT_GENERATE_FOR_TYPE \n"
#endif // USE_ORIGINAL_MACROS_SHOULD_BE_FALSE
