#include "gpuRandom.hpp"


#define DEBUG



// copy output = t(input)
// I don't think it will work with float.
//[[Rcpp::export]]
Rcpp::NumericMatrix copyToCpu(
    Rcpp::S4 input
  ) {
  
  Rcpp::traits::input_parameter< std::string >::type classVarR(RCPP_GET_CLASS(input));
  std::string precision_type = (std::string) classVarR;
  
  if(precision_type != "dvclMatrix") {
        Rcpp::warning("class of input must be dvclMatrix");
  }
  
  // TO DO: use async_copy to copy by column and vector_iterator

  const bool BisVCL=1;
  const int ctx_id = INTEGER(input.slot(".context_index"))[0]-1;
  std::shared_ptr<viennacl::matrix_base<double> > xMat = getVCLptr<double>(input.slot("address"), BisVCL, ctx_id);
  viennacl::vector_base<double> xVec((*xMat).handle(), (*xMat).internal_size(), 0, 1);
  
  const int Nrow = (*xMat).size1(), Ncol = (*xMat).size2(), NcolInt = (*xMat).internal_size2();


//  viennacl::vector_base<double> xRow;
Rcpp::NumericMatrix output(Nrow, Ncol);  
Rcpp::NumericVector outputVec(output), oneRow(Ncol);
double *outputColP = &(outputVec[0]), *oneRowP = &(oneRow[0]);


  int Drow;

//      viennacl::vector_iterator<double,1> xIter = xVec.begin();
      
      
    int Nbytes =sizeof(double)*Ncol, NbytesPerRow = sizeof(double)*NcolInt;
      
    viennacl::backend::mem_handle xHandle =(*xMat).handle(); 
    viennacl::ocl::handle< cl_mem > xHandleOcl = xHandle.opencl_handle();
    viennacl::ocl::context & memory_context = const_cast<viennacl::ocl::context &>(xHandleOcl.context());
 //   cl_command_queue  theQueue = memory_context.get_queue().handle().get();
//    cl_mem theMem = xHandleOcl.get();
#ifdef DEBUG
    Rcpp::Rcout << "start" << "\n";
#endif    
    
    for(Drow = 0; Drow < Nrow; Drow++) {

//      viennacl::backend::memory_read( (*xMat).handle(),
      viennacl::backend::opencl::memory_read( xHandleOcl,
                                      Drow*NbytesPerRow, 
                                      Nbytes, 
                                     oneRowP,
                                     0 //async
                             );
      output(Drow, Rcpp::_) = oneRow;
      
#ifdef DEBUG
      if(Drow < 10 ) {
        Rcpp::Rcout << Drow << " ";
      }
#endif    
      
#ifdef UNDEF            
        viennacl::vector_iterator<double,1> xIterHere = xIter + (Drow*NcolInt);
            viennacl::vector_iterator<double,1> xIterHereEnd = xIterHere + Ncol;
      viennacl::fast_copy(
      xIterHere, 
      xIterHereEnd,
      outputColP + Drow*Ncol
    );
#endif
    }
    Rcpp::Rcout << "done" << "\n";
    
      
      return(output); 
//  return(output); 
}