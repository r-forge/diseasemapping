


template <typename T> 
std::string multiplyLowerTypeString() {

  std::string typeString = openclTypeString<T>();
  std::string result = "";

  if(typeString == "double") {
  	result += "\n#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n\n";
  }
  result = result + 
"__kernel void cholGpu("
"	__global double *A," 
"	__global double *diag,"
"	__local double *diagLocal,"
"	const int N, "
"	const int Npad,"
"	const int NlocalStorage){"
// first dimension is rows, second dimension is work items
// doing the same row
"	const int Dglobal0 = get_global_id(0);"
"	const int Dglobal1 = get_global_id(1);"
"	const int Nsize0 = get_global_size(0);"
"	const int Nsize1 = get_global_size(1);"
"	const int NlocalSize = get_local_size(1);"
 
"	int Dcol=0, Dcolm1;"
"	int DcolNpad = Dcol*Npad;"
"	int maxDcolNlocalStorage = min(Dcol, NlocalStorage);"

"	__private int Drow, Dk, DrowDk;"
"	__private double DL, Ddiag;\n"

"for(Dcol = 0; Dcol < N; Dcol++){\n"
"   Dcolm1 = Dcol-1;\n"
"   DcolNpad = Dcol*Npad;\n"
// diagonal entry, use only first group
" if(Dglobal0 = 0){\n"
"	Ddiag = 0.0;\n"
"	for(Dk = Dglobal1; Dk < Dcol; Dk += Nsize1) {\n"
"		DL = A[Dcol + Dk * Npad];"
"		DL *= DL;"
"		diagLocal[Dlocal] += diag[Dk] * DL;"
"	}\n"
"	}\n" // Dglobal0 == 0
"	barrier(CLK_LOCAL_MEM_FENCE);\n"
// add 'em up to get diagonal
"	if(Dglobal0==0 & Dglobal1==0) {\n"
"		Ddiag = 0.0;\n"
"		for(Dk = 0; Dk < NlocalSize; Dk++) {\n"
"			Ddiag += diagLocal[Dk];\n"
"		}\n"
"		diag[Dcol] = A[Dcol + DcolNpad] - Ddiag;\n"
"	}\n" // Dglobal1 == 0
"	barrier(CLK_LOCAL_MEM_FENCE);\n"
 
	// diagLocal =  diag * A[Dcol, ]
	// TO DO: use less memory for diagLocal
"	for(Dk = Dlocal; Dk < maxDcolNlocalStorage; Dk += NlocalSize) {\n"
"		DL = A[Dcol + Dk * Npad];\n"
"		diagLocal[Dk] = diag[Dk] * DL;\n"
"	}\n"
"	barrier(CLK_LOCAL_MEM_FENCE);\n" 

"	for(Drow = (N-Dglobal-1); Drow > Dcol; Drow -= Nsize) {\n"

"		DL = A[Drow + DcolNpad];\n"

"		DrowDk=Dcolm1*Npad;\n"
"		for(Dk = Dcolm1; Dk >= 0; Dk--) {\n"
//			DL -= A[DrowDk] * diagLocal[Dk];"
"			DrowDk -= Npad;"
"			DL -= A[Drow + Dk * Npad] * diagLocal[Dk];"
			// DL -= A[Drow + Dk * Npad] * A[Dcol + DkNpad] * diag[Dk];"
"		}\n" // Dk
"		A[Drow + DcolNpad] = DL / Ddiag;"

"	}\n" // Drow
	// Ddiag is now diag[Dcol]
	// copy it to global memory
"	if(Dglobal == 0) {"
"		diag[Dcol] = Ddiag;"
"		A[Dcol + DcolNpad] = 1.0;"
"	}\n"
"}\n" // Dcol loop
"}";
