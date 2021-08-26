#include "gpuRandom.hpp"


//#define DEBUGKERNEL

// TO DO: use async_copy to copy by column and vector_iterator
template<typename T> 
void copyToCpuTemplated(
    Rcpp::NumericMatrix output,
    Rcpp::S4 input) {
  
  const bool BisVCL=1;
  const int ctx_id = INTEGER(input.slot(".context_index"))[0]-1;
  std::shared_ptr<viennacl::matrix_base<T> > xMat = getVCLptr<T>(input.slot("address"), BisVCL, ctx_id);

  viennacl::vector_base<T> x((*xMat).handle(), (*xMat).internal_size(), 0, 1);
  
  
  int Drow=0;
  const int Nrow = (*xMat).size1(), Ncol = (*xMat).size2(), NcolInt = (*xMat).internal_size2();
  for(Drow = 0; Drow < Nrow; Drow++) {
    viennacl::fast_copy(
      x.begin() + Drow*NcolInt,// + Drow*NcolInt, 
      x.begin() + Drow*NcolInt + Ncol,// + Drow*NcolInt+Ncol, 
      output.begin() +Drow*Ncol
    );
  }
  
}

// copy output = t(input)
// I don't think it will work with float.
//[[Rcpp::export]]
void copyToCpu(
    Rcpp::NumericMatrix output,
    Rcpp::S4 input
  ) {
  
  
  Rcpp::traits::input_parameter< std::string >::type classVarR(RCPP_GET_CLASS(input));
  std::string precision_type = (std::string) classVarR;
  
  if(precision_type == "fvclMatrix") {
    copyToCpuTemplated<float>(
      output, input);
  } else if (precision_type == "dvclMatrix") {
    copyToCpuTemplated<double>(
      output, input);
  } else {
    Rcpp::warning("class of input must be fvclMatrix or dvclMatrix");
  }
}