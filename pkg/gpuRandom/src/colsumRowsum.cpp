

template <typename T> 
std::string colsumRowsumString(
    const int Nrow, 
    const int Ncol,
    const int NpadCol) { 
  
  std::string typeString = openclTypeString<T>();
  std::string result = "";
  
  if(typeString == "double") {
    
    result += "\n#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n"
    
  }

  
  result +=
    "\n#define Nrow " + std::to_string(Nrow) + "\n"    
    "#define Ncol " + std::to_string(Ncol) + "\n"
    "#define NpadCol " + std::to_string(NpadCol) + "\n";    
  
    result += mrg31k3pString();
  
  result += 
    "\n\n__kernel void colsumRowsum(\n"
    "  __global " + typeString + "* x,"  
    "  __global " + typeString + "* rowSum,"  
    "  __global " + typeString + "* colSum){\n\n";  
  
  
 
  result += "int Drow, Dcol, Dindex;\n";
  result += typeString + " result;\n";
  
  result + "if(get_global_id(1) { // sum columns\n"
  "  for(Drow = get_global_id(0); Drow < Nrow; Drow++{\n"
  "    result = 0.0;\n"
  "    for(Dcol = 0, Dindex = Drow*NpadCol; Dcol < Ncol; Dcol++, Dindex++){\n"
  "       result += x[Dindex];\n"    
  "    } // end loop through columns\n"
  "    colSum[Dcol] = result;\n"
  "  } // end loop through rows\n"
  "} else { // sum rows\n"
  "  for(Dcol = get_global_id(0); Dcol < Ncol; Dcol++){\n"
  "    result = 0.0;\n"
  "    for(Drow = 0, Dindex = Dcol; Drow < Nrow; Drow++, Dindex+= NpadCol){\n"
  "       result += x[Dindex];\n"    
  "    } // end loop through columns\n"
  "    rowSum[Drow] = result;\n"
  "  } // end loop through rows\n"

  "}\n\n";
  result += 
    "}//kernel\n";
  
  return(result);
}

