
#include "viennacl/vector_proxy.hpp"
#include "viennacl/matrix_proxy.hpp"


void maternBatchBackend(Rcpp::S4 var,
    Rcpp::S4 coords,
    Rcpp::S4 param, 
    Rcpp::IntegerVector Nglobal,
    Rcpp::IntegerVector Nlocal,
    int startrow,   
    int numberofrows,
    int verbose=0) 
    
    
void cholBatchBackend(Rcpp::S4 A,
        Rcpp::S4 D,
        Rcpp::IntegerVector Astartend,
        Rcpp::IntegerVector Dstartend,
        const int numbatchD,
        std::vector<int> Nglobal,
        std::vector<int> Nlocal,
        std::vector<int> NlocalCache)     
    
    
void rowsumBackend(
        Rcpp::S4  xR, 
        Rcpp::S4  SumR,
        std::string type,
        int log)
        
        
SEXP backsolveBatchBackend(
            Rcpp::S4 C,
            Rcpp::S4 A,
            Rcpp::S4 B,
            Rcpp::IntegerVector Cstartend,
            Rcpp::IntegerVector Astartend,
            Rcpp::IntegerVector Bstartend,
            const int numbatchB,
            const int diagIsOne,
            Rcpp::IntegerVector Nglobal,
            Rcpp::IntegerVector Nlocal, 
            const int NlocalCache)       
        
        
SEXP crossprodBatchBackend(
            Rcpp::S4 C,
            Rcpp::S4 A,
            Rcpp::S4 D,
            const int invertD,
            Rcpp::IntegerVector Cstartend,
            Rcpp::IntegerVector Astartend,
            Rcpp::IntegerVector Dstartend, 
            Rcpp::IntegerVector Nglobal,
            Rcpp::IntegerVector Nlocal, 
            const int NlocalCache)        
      
        
SEXP gemmBatch2backend(Rcpp::S4 A,
                Rcpp::S4 B,  
                Rcpp::S4 C,
                Rcpp::IntegerVector transposeABC,  
                Rcpp::IntegerVector submatrixA,
                Rcpp::IntegerVector submatrixB,
                Rcpp::IntegerVector submatrixC, 
                Rcpp::IntegerVector batches, 
                Rcpp::IntegerVector workgroupSize,   
                Rcpp::IntegerVector NlocalCache,
                const int verbose)        
        
        
void matrix_vector_sumBackend(Rcpp::S4 matrixR,
                Rcpp::S4 vectorR,
                Rcpp::S4 sumR,
                const int byrow,
                Rcpp::IntegerVector numWorkItems)  
                
                
                
                
                
                
                
                
                
                
                
        