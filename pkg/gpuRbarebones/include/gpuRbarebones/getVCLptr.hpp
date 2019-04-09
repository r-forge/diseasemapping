
#ifndef GET_VCL_PTR_HPP
#define GET_VCL_PTR_HPP

#include "gpuRbarebones/dynVCLMat.hpp"
#include "gpuRbarebones/dynVCLVec.hpp"


template<typename T>
std::shared_ptr<viennacl::matrix<T> >
getVCLptr(
    SEXP ptr_,
    const bool isVCL,
    const int ctx_id
){
    std::shared_ptr<viennacl::matrix<T> > vclptr;
    
    if(isVCL){
        Rcpp::XPtr<dynVCLMat<T> > ptr(ptr_);
        vclptr = ptr->sharedPtr();
    }

    
    return vclptr;
}

template<typename T>
std::shared_ptr<viennacl::vector_base<T> >
getVCLVecptr(
    SEXP ptr_,
    const bool isVCL,
    const int ctx_id
){
    std::shared_ptr<viennacl::vector_base<T> > vclptr;
    
    if(isVCL){
        Rcpp::XPtr<dynVCLVec<T> > ptr(ptr_);
        vclptr = ptr->sharedPtr();
    }
    
    
    return vclptr;
}

#endif
