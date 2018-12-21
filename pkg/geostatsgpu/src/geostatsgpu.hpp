
//#include <string>
#include <Rcpp.h>
//#include "gpuR/getVCLptr.hpp"

//#include "dynVCLMatGeostatsgpu.hpp"
#include "gpuR/dynVCLMat.hpp"
#include "dynVCLVecGeostatsgpu.hpp"

#include "viennacl/ocl/backend.hpp"
#include "viennacl/linalg/sum.hpp"
#include "viennacl/linalg/lu.hpp"

#define GSL_DBL_EPSILON 2.2204460492503131e-16
#define GSL_SQRT_DBL_MAX 1.3407807929942596e+154
extern "C" void Rtemme_gamma(double *nu, double * g_1pnu, double * g_1mnu, double *g1, double *g2);



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
