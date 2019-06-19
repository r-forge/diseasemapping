template <typename T> 
int sizeOfReal() {
  return(-1);}

template <> int sizeOfReal<double>(){
  return(sizeof(cl_double));}

template <> int sizeOfReal<float>(){
  return(sizeof(cl_float));}



template <typename T> 
std::string openclTypeString() {
  return("undefined");}

template <> std::string openclTypeString<double>(){
  std::string result = "double";
  return(result);
}

template <> std::string openclTypeString<float>(){
  std::string result = "float";
  return(result);
}

template <> std::string openclTypeString<int>(){
  std::string result = "int";
  return(result);
}
