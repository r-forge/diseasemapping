
#include <Rcpp.h>
#include <string>

#include <dynMatrix/dynVCLMatGeostatsgpu.hpp>
#include <dynMatrix/dynVCLVecGeostatsgpu.hpp>
#include "viennacl/linalg/sum.hpp"
#include "viennacl/ocl/backend.hpp"


template <typename T> std::string maternBatchKernelString();
std::string cholBatchKernelString();
std::string backsolveBatchString();
std::string crossprodBatchString();
std::string gemmBatch2String();
std::string matrix_plus_vectorString();
std::string matrix_plus_scalarString();


//////////////////////////////////////////////////////////////////////////////////////////////////////
//#include "viennacl/linalg/lu.hpp"

extern "C" void Rtemme_gamma(double *nu, double * g_1pnu, double * g_1mnu, double *g1, double *g2);

//template<typename T> double luT(viennacl::matrix<T> &vclX, viennacl::vector_base<T> &vclD);

template <typename T> T maternClEpsilon();




#include "viennacl/vector_proxy.hpp"
#include "viennacl/matrix_proxy.hpp"





template<typename T> 
void maternBatchVcl(
    viennacl::matrix<T> &vclVar, // Nmat columns N^2 rows
    viennacl::matrix<T> &vclCoords, // 2 columns
    viennacl::matrix<T> &param, // Nmat rows, 22 columns
    Rcpp::IntegerVector numWorkItems,
    Rcpp::IntegerVector numLocalItems,	
    const int ctx_id,
    int startrow,   // new added
    int numberofrows,
    int verbose);
    
    
    template <typename T> 
    int cholBatchVcl(
        viennacl::matrix<T> &A,
        viennacl::matrix<T> &D,
        Rcpp::IntegerVector Astartend,
        Rcpp::IntegerVector Dstartend,  
        const int numbatchD,
        Rcpp::IntegerVector Nglobal,
        Rcpp::IntegerVector Nlocal, 
        Rcpp::IntegerVector NlocalCache,
        const int ctx_id);
    
    
    template <typename T>
    void rowsum(
        viennacl::matrix<T> &x,
        viennacl::vector_base<T> &Sum,
        std::string type,
        int log);
        
        
        template <typename T> 
        void backsolveBatch(
            viennacl::matrix<T> &C,
            viennacl::matrix<T> &A,  //must be batches of square matrices
            viennacl::matrix<T> &B,
            Rcpp::IntegerVector Cstartend,
            Rcpp::IntegerVector Astartend, //square matrices
            Rcpp::IntegerVector Bstartend, 
            const int numbatchB,
            const int diagIsOne,
            Rcpp::IntegerVector Nglobal,
            Rcpp::IntegerVector Nlocal,
            const int NlocalCache,
            const int ctx_id);       
        
        
        template <typename T> 
        void crossprodBatch(
            viennacl::matrix<T> &C,  // must be a batch of square matrices 
            viennacl::matrix<T> &A,
            viennacl::matrix<T> &D,
            const int invertD,
            Rcpp::IntegerVector Cstartend,
            Rcpp::IntegerVector Astartend,
            Rcpp::IntegerVector Dstartend,  
            Rcpp::IntegerVector Nglobal,
            Rcpp::IntegerVector Nlocal,
            const int NlocalCache, 
            const int ctx_id);      
      
        
        template <typename T> 
        int gemmBatch2(
            viennacl::matrix<T> &A,
            viennacl::matrix<T> &B,
            viennacl::matrix<T> &C,
            Rcpp::IntegerVector transposeABC,  
            Rcpp::IntegerVector submatrixA,
            Rcpp::IntegerVector submatrixB,
            Rcpp::IntegerVector submatrixC,  
            Rcpp::IntegerVector batches, 
            Rcpp::IntegerVector workgroupSize,
            Rcpp::IntegerVector NlocalCache, 
            const int verbose,
            const int ctx_id); 
        
        
        
template<typename T> 
void matrix_vector_sum(
     viennacl::matrix<T> &matrix,// viennacl::vector_base<int>  rowSum, viennacl::vector_base<int>  colSum,  
     viennacl::vector_base<T> &vector,
     viennacl::matrix<T> &sum,
     const int byrow,
     Rcpp::IntegerVector numWorkItems,
     int ctx_id);            
                
                
template<typename T> 
void matrix_scalar_sum(
      viennacl::matrix<T> &matrix,// viennacl::vector_base<int>  rowSum, viennacl::vector_base<int>  colSum,  
      T  value,         //viennacl::scalar<T> value,     
      viennacl::matrix<T> &sum,
      Rcpp::IntegerVector numWorkItems,
      int ctx_id);                
                
                
                
                
                
        