
#include "lgmlikFit.hpp"
//#include "gpuRandom.hpp"
//#include<Rmath.h>

using namespace Rcpp;
using namespace viennacl; 
using namespace viennacl::linalg;


//#define M_PI  3.141592653589793238463
/*
 const double PI  =3.141592653589793238463;
 const float  PI_F=3.14159265358979f;*/
/*#define M_LN_SQRT_2PI 0.918938533204672741780329736406  /* log(sqrt(2*pi))
 == log(2*pi)/2 */


template<typename T> 
void likfitGpu_0(viennacl::matrix<T> &Vbatch,
                 viennacl::matrix<T> &coordsGpu, 
                 viennacl::matrix_range<viennacl::matrix<T> > &paramsBatch, 
                 viennacl::matrix<T> &diagMat, 
                 viennacl::vector_range<viennacl::vector<T> > &logD,
                 viennacl::matrix<T> &ab,
                 viennacl::matrix<T> &yX, //y1,y2,y3,X            
                 viennacl::matrix<T> &temp2,
                 viennacl::matrix_range<viennacl::matrix<T> > &aTDa,   // ssqY
                 viennacl::matrix<T> &betas, //a vclmatrix  //given by the user or provided from formula, default=null
                 viennacl::matrix<T> &temp00,
                 viennacl::matrix<T> &temp0,
                 viennacl::matrix<T> &temp1,
                 viennacl::matrix<T> &ssqBeta0,
                 viennacl::matrix_range<viennacl::matrix<T> > &ssqBeta, //  
                 viennacl::matrix<T> &one,          
                 viennacl::matrix<T> &diagP,
                 viennacl::matrix<T> &temp3,
                 viennacl::matrix<T> &nine0,
                 viennacl::matrix_range<viennacl::matrix<T> > &nine,      // ssqX
                 viennacl::matrix<T> &two,
                 viennacl::vector_range<viennacl::vector<T> > &logP,  // matrix logPlogD
                 viennacl::matrix<T> &Qinverse,
                 viennacl::matrix<T> &identity,
                 viennacl::matrix<T> &QPQinverse,
                 viennacl::matrix_range<viennacl::matrix<T>> &betahat,  // p*rowbatch   colbatch                            
                 viennacl::matrix_range<viennacl::matrix<T>> &variances,   // must be rowbatch * colbatch matrix !!!
                 viennacl::vector_base<T> &logD_plusP,
                 viennacl::matrix<T> &form_temp0,
                 viennacl::matrix<T> &form_temp,
                 viennacl::matrix<T> &form_temp1,
                 viennacl::matrix<T> &form_temp2,
                 viennacl::matrix<T> &jacobian,    
                 viennacl::matrix_range<viennacl::matrix<T> > &LogLik,  
                 int n, 
                 int p, 
                 int rowbatch, 
                 int colbatch,   
                 int form, // c(loglik=1, ml=2, mlFixSigma=3, mlFixBeta=4, reml=5, remlPro=6)
                 Rcpp::IntegerVector workgroupSize,
                 Rcpp::IntegerVector localSize,
                 Rcpp::IntegerVector localSizechol,
                 Rcpp::IntegerVector NlocalCache,
                 const int ctx_id){
  
  const T nlog2pi  = (T) n*2*M_LN_SQRT_2PI;   //  M_LN_SQRT_2PI == log(2*pi)/2  
  // logP = column 1 of logPlogD, logD = column2 of logPlogD
  
  
  //if((form == 1 && betas.size1==betas.size2==0) || (form == 3 && betas.size1==betas.size2==0) )  
  
  //  Rcpp::Rcout << "\n" << "need betas" << "\n";
  //get matern matrix
  maternBatchVcl(Vbatch, coordsGpu, paramsBatch, workgroupSize, localSize, ctx_id, 0, rowbatch, 0);
  //#Vbatch=LDL^T, cholesky decomposition
  Rcpp::IntegerVector Astartend =IntegerVector::create(0,n,0,n);
  Rcpp::IntegerVector Dstartend =IntegerVector::create(0,1,0,n);
  cholBatchVcl(Vbatch, diagMat, Astartend, Dstartend, rowbatch, workgroupSize, localSizechol, NlocalCache, ctx_id);
  
  
  //logD <- apply(log(diagMat),1,sum)
  // half log determinant of V
  rowsum(diagMat, logD, "row", 1);
  //variances <- vclMatrix(matrix(paramsBatch[,3], nrow=rowbatch, ncol=colbatch, byrow=FALSE), type=type)     # this has to be a vclMatrix not vector
  
  
  //#L(a1,a2,a3, b) = (y1,y2,y3, X)
  Rcpp::IntegerVector Cstartend =IntegerVector::create(0, n,0, colbatch+p);
  //Astartend = IntegerVector::create(0, n,0, n);
  Rcpp::IntegerVector Bstartend =IntegerVector::create(0,n,0,colbatch+p);  
  backsolveBatch(ab, Vbatch, yX, Cstartend, Astartend, Bstartend, 1, 1, workgroupSize, localSize, NlocalCache[0], ctx_id);
  
  
  
  //# temp2 = (ab)^T * D^(-1) *ab  = (Y X)^T V^{-1} (YX)
  Cstartend=IntegerVector::create(0,colbatch+p,0,colbatch+p);
  Astartend=IntegerVector::create(0,n,0,colbatch+p);
  //Dstartend=IntegerVector::create(0,1,0,n);
  crossprodBatch(temp2, ab, diagMat, 1, Cstartend, Astartend, Dstartend,  workgroupSize, localSize, NlocalCache[0], ctx_id);
  
  /*extract a^TDa from temp2       1x3
   for (j in 1:colbatch){
   for (i in 1:rowbatch){
   aTDa[i,j]= temp2[(i-1)*ncol(ab)+j, j]
   }
   }*/
  
  /*
   
   SAVE SSQY
   
   */
  // extract a^TDa from temp2       1x3
  // multiple datasets, save as aTDa or ssqY
  viennacl::slice s2(0, 1, colbatch);   // {0,1,2}  colbatch=3    ？？？？？？？？？？？
  for (int i = 0; i < rowbatch; i++){// square sub_matrices
    viennacl::slice s1((colbatch+p)*i, 1, colbatch); 
    matrix_slice<viennacl::matrix<T> >  temp2_sub(temp2, s1, s2);  //project(temp2, r, r); //returns a matrix_range as above
    viennacl::vector<T> diags = viennacl::row(aTDa, i);
    diags = viennacl::diag(temp2_sub);
  }
  
  /*
   cholesky of X^T V^(-1) X
   */
  // to get hat_sigma^2 when beta is not given by users
  // b^T * D^(-1) * b = Q * P * Q^T, cholesky of a subset (right bottom) of temp2
  Astartend =IntegerVector::create(colbatch, p, colbatch, p);
  Dstartend =IntegerVector::create(0, 1, 0, p);
  cholBatchVcl( temp2, diagP, Astartend, Dstartend,  rowbatch, workgroupSize, localSizechol, NlocalCache, ctx_id);
  // store determinant for REML  
  rowsum(diagP, logP, "row",1L);  
  
  /* 
   SSQX
   */
  // ssqx = (X betahat)^T V^(-1) X betahat = nine
  // temp3 = Q^(-1) * (b^T * D^(-1) *a), backsolve for temp3    2 by 1
  Cstartend =IntegerVector::create(0,p,0,colbatch);  //Astartend = IntegerVector::create(colbatch, p, colbatch, p);
  Bstartend =IntegerVector::create(colbatch, p, 0, colbatch);  
  backsolveBatch(temp3, temp2, temp2, Cstartend, Astartend, Bstartend, rowbatch, 1, workgroupSize, localSize, NlocalCache[0], ctx_id);
  
  
  
  // nine0 = temp3^T * P^(-1) * temp3,  four 1 by 1 matrices      3x3
  Cstartend = IntegerVector::create(0, colbatch, 0, colbatch);
  Astartend = IntegerVector::create(0, p, 0, colbatch); //Dstartend =IntegerVector::create(0, 1, 0, p);  
  crossprodBatch(nine0, temp3, diagP, 1,  Cstartend, Astartend, Dstartend,  workgroupSize, localSize, NlocalCache[0], ctx_id);
  // TO DO (long list), save SSQX matrix, row batches and column batches
  
  // save diagonals of ssqX, for standard errors of beta hat
  // extract needed cells from nine0       1x3
  for (int i = 0; i < rowbatch; i++){
    viennacl::slice s1(colbatch*i, 1, colbatch); 
    matrix_slice<viennacl::matrix<T> >  nine0_sub(nine0, s1, s2);  //project(temp2, r, r); //returns a matrix_range as above
    viennacl::vector<T> diag4 = viennacl::row(nine, i);
    diag4 = viennacl::diag(nine0_sub);
  }
  
  
  // betahat = (X^T V^(-1) X)^(-1)  X^T V^(-1) Y
  // Q^T P Q = (X^T V^(-1) X)
  // calculate betahat = Q^(-T) p^(-1) Q^(-1) * b^T D^(-1)a
  // first: get Q^(-1)    pxp
  Cstartend = IntegerVector::create(0, p, 0, p);
  Astartend = IntegerVector::create(colbatch, p, colbatch, p);
  Bstartend = IntegerVector::create(0, p, 0, p);    
  backsolveBatch(Qinverse, temp2, identity,Cstartend, Astartend, Bstartend, 1,  1, workgroupSize, localSize, NlocalCache[0], ctx_id);
  
  
  // Q^(-T) P^(-1) Q^(-1) = QPQinverse = SSQX^(-1)
  //Cstartend = IntegerVector::create(0, p, 0, p);
  Astartend = IntegerVector::create(0, p, 0, p);
  //Dstartend = IntegerVector::create(0, 1, 0, p);    
  crossprodBatch(QPQinverse, Qinverse, diagP, 1,Cstartend, Astartend, Dstartend,  workgroupSize, localSize, NlocalCache[0], ctx_id);
  
  
  // QPQinverse * b^T D^(-1)a = betahat
  //Rcpp::IntegerVector transposeABC = IntegerVector::create(0,0,0);
  Rcpp::IntegerVector transposeABC = IntegerVector::create(0,0,0);
  Rcpp::IntegerVector batches = IntegerVector::create(rowbatch,1,0,0,0,0); //batches: nRow, nCol, recycleArow, recycleAcol, recycleBrow, col
  Rcpp::IntegerVector workgroupSize_gemm = IntegerVector::create(workgroupSize[0],workgroupSize[1],workgroupSize[2],localSize[0],localSize[1],localSize[2]);
  Rcpp::IntegerVector NlocalCache_gemm = IntegerVector::create(NlocalCache[0],NlocalCache[0]);
  
  Astartend=IntegerVector::create(0,p,p,0,p,p);
  Bstartend=IntegerVector::create(colbatch, p,(colbatch+p), 0, colbatch, (colbatch+p));
  Cstartend=IntegerVector::create(0,p,p,0,colbatch,colbatch);
  gemmBatch2(QPQinverse, temp2, betahat,transposeABC,  Astartend, Bstartend, Cstartend,  batches, workgroupSize_gemm, NlocalCache_gemm, 0, ctx_id);
  
  
  
  
  // resid^T V^(-1) resid, resid = Y - X betahat 
  // nine is ssqx
  two = aTDa - nine;
  // log ssqResid
  form_temp0 = element_log(two);  // overwritten later, keep here for debugging
  LogLik = form_temp0;
  
  /* important things saved are
   logLik is logSsqResid
   aTDa = SSQY 
   nine = SSQX
   logD, logP, determinants, return, matrix with Nparam rows, two columns
   */
  
#ifdef UNDEF
  
  /*
   // Extract 5-th row of A, then overwrite with 6-th diagonal:
   viennacl::vector<T> r = viennacl::row(A, 4);
   r = viennacl::row(A, 5);
   viennacl::vector<T> temp2_diag = viennacl::diag(temp2);
   viennacl::range r(0, colbatch); 
   vector_range<vector<T> > aTDa(temp2_diag, r); */
  
  
  batches = IntegerVector::create(rowbatch,1,0,0,1,0); //batches: nRow, nCol, recycleArow, recycleAcol, recycleBrow, col
  
  if(form == 1 || form == 3){      // beta is supplied 
    // a^TD^(-1)b * beta = temp00
    Astartend=IntegerVector::create(0,colbatch,(colbatch+p),colbatch,p,(colbatch+p));
    Bstartend=IntegerVector::create(0,p,p,0,colbatch,colbatch);
    Cstartend=IntegerVector::create(0,colbatch,colbatch,0,colbatch,colbatch);
    gemmBatch2(temp2, betas, temp00, transposeABC,  Astartend, Bstartend, Cstartend, batches, workgroupSize_gemm, NlocalCache_gemm, 0, ctx_id);
    
    
    // extract a^TD^(-1)b * beta from temp00      1x3
    for (int i = 0; i < rowbatch; i++){
      viennacl::slice s1(colbatch*i, 1, colbatch); 
      matrix_slice<viennacl::matrix<T> >  temp00_sub(temp00, s1, s2);  //project(temp2, r, r); //returns a matrix_range as above
      viennacl::vector<T> diag2 = viennacl::row(temp0, i);
      diag2 = viennacl::diag(temp00_sub);
    }
    
    
    
    // b * beta = temp1 (n x colbatch), to get ssqBeta
    Astartend=IntegerVector::create(0,n,n,colbatch,p,(colbatch+p));
    //Bstartend=IntegerVector::create(0,p,p,0,colbatch,colbatch);
    Cstartend=IntegerVector::create(0,n,n,0,colbatch,colbatch);
    //Rcpp::IntegerVector batches = IntegerVector::create(rowbatch,1,0,0,1,0); 
    gemmBatch2( ab, betas, temp1, transposeABC,  Astartend, Bstartend, Cstartend, batches,  workgroupSize_gemm, NlocalCache_gemm, 0, ctx_id);
    
    
    
    
    // ssqBeta0 = temp1^T D^(-1) temp1 = beta T*bT D^(-1) b*beta    3x3
    Cstartend=IntegerVector::create(0, colbatch, 0, colbatch);
    Astartend=IntegerVector::create(0, n, 0, colbatch);
    //Dstartend=IntegerVector::create(0, 1, 0, n);
    crossprodBatch( ssqBeta0, temp1, diagMat,  TRUE,  Cstartend, Astartend, Dstartend,  workgroupSize, localSize, NlocalCache[0], ctx_id);
    
    
    
    // extract beta T*bT D^(-1) b*beta from ssqBeta0       1x3
    for (int i = 0; i < rowbatch; i++){
      viennacl::slice s1(colbatch*i, 1, colbatch); 
      matrix_slice<viennacl::matrix<T> >  ssqBeta0_sub(ssqBeta0, s1, s2);  //project(temp2, r, r); //returns a matrix_range as above
      viennacl::vector<T> diag3 = viennacl::row(ssqBeta, i);
      diag3 = viennacl::diag(ssqBeta0_sub);
    }
    
    // to get one, see the notes 
    one = aTDa - 2*temp0 + ssqBeta;      // 1x3
  }
  
  
  
  
  
  
  
  if(form == 1 || form == 4){     // variances are provided
    // get form_temp <- n*log(variances) + logD    , form_temp = n*log(sigma^2)+log |D|
    viennacl::matrix<T> nlog_variances(rowbatch, colbatch);  
    nlog_variances = n * element_log(variances);
    matrix_vector_sum( nlog_variances, logD, form_temp, 1, workgroupSize, ctx_id);  // form_temp
    
  } 
  
  if (form ==5 || form ==6){
    logD_plusP = logD + logP;   
  }
  
  
  
  if (form==1){ //loglik
    //form_temp + one/variances + n*log(2*PI)+jacobian 
    //form_temp0 = form_temp + element_div(one, variances);
    //matrix_vector_sum(form_temp0, jacobian, form_temp2, 0L, workgroupSize, ctx_id);
    form_temp0 = form_temp + element_div(one, variances) + jacobian;
    
    matrix_scalar_sum(form_temp0, nlog2pi, form_temp1, workgroupSize,ctx_id);
    
    LogLik = -0.5*form_temp1;
    
  }else if(form ==2){ //ml  =n*log(two/n) +logD +jacobian + n*log(2*PI) + n 
    form_temp0 = n * element_log(two/n);
    
    matrix_vector_sum(form_temp0, logD, form_temp, 1L, workgroupSize, ctx_id);   
    //matrix_vector_sum(form_temp, jacobian, form_temp2, 0L, workgroupSize, ctx_id);
    
    form_temp2 = form_temp + jacobian;
    
    matrix_scalar_sum(form_temp2, nlog2pi+n, form_temp1, workgroupSize,ctx_id);
    
    LogLik = -0.5*form_temp1;
    
  }else if(form==3){ // mlFixSigma
    //n*log(one/n)+logD + n*log(2*PI) + n+jacobian
    form_temp0 = n * element_log(one/n);
    
    matrix_vector_sum(form_temp0, logD, form_temp, 1L, workgroupSize, ctx_id);
    //matrix_vector_sum(form_temp, jacobian, form_temp2, 0L, workgroupSize, ctx_id);
    form_temp2 = form_temp + jacobian;
    
    matrix_scalar_sum(form_temp2, nlog2pi+ n, form_temp1, workgroupSize,ctx_id);
    
    LogLik = -0.5*form_temp1;
    
  }else if(form==4){ //mlFixBeta
    //form_temp + two/variances + n*log(2*PI)+jacobian 
    //form_temp0 = form_temp + element_div(two, variances);
    //matrix_vector_sum(form_temp0, jacobian, form_temp2, 0L, workgroupSize, ctx_id);
    form_temp0 = form_temp + element_div(two, variances) + jacobian;
    
    matrix_scalar_sum(form_temp0, nlog2pi, form_temp1, workgroupSize,ctx_id);
    
    LogLik = -0.5*form_temp1;
    
  }else if(form==5){ //reml
    //form_temp = (n-p)*log(variances) + logD + logP
    //minusTwoLogLik=form_temp + two/variances + n*log(2*PI)+jacobian
    form_temp0 = (n-p) * element_log(variances);
    matrix_vector_sum(form_temp0, logD_plusP, form_temp, 1L, workgroupSize, ctx_id);    
    //viennacl::matrix<T> plus_two_div_variances = form_temp + element_div(two, variances);
    //matrix_vector_sum(plus_two_div_variances, jacobian, form_temp2, 0L, workgroupSize, ctx_id);
    form_temp2 = form_temp + element_div(two, variances) + jacobian;
    
    matrix_scalar_sum(form_temp2, nlog2pi, form_temp1, workgroupSize, ctx_id);
    
    
    LogLik = -0.5 * form_temp1;
    
  }else if(form==6){ // remlPro
    // minusTwoLogLik= (n-p)*log(two/(n-p)) + logD + logP + jacobian + n*log(2*PI) + n-p , 
    form_temp0 = (n-p) * element_log(two/(n-p));
    matrix_vector_sum(form_temp0, logD_plusP, form_temp, 1L, workgroupSize, ctx_id);   
    //matrix_vector_sum(form_temp, jacobian, form_temp2, 0L, workgroupSize, ctx_id);
    form_temp2 =form_temp+ jacobian;
    
    matrix_scalar_sum(form_temp2, nlog2pi + n-p, form_temp1, workgroupSize, ctx_id);
    
    LogLik = -0.5*form_temp1;
  } 
  
#endif 
  
}










template<typename T> 
void likfitGpu_1(viennacl::matrix<T> &coordsGpu, 
                 viennacl::matrix<T> &bigparamsBatch, 
                 viennacl::matrix<T> &yX, //y1,y2,y0,X,   already transformed            
                 viennacl::matrix<T> &betas, //a vclmatrix  //given by the user or provided from formula, default=null                          
                 viennacl::matrix<T> &bigvariances,   // must be rowbatch * colbatch matrix !!!
                 viennacl::matrix<T> &jacobian,    // (groupsize, colbatch)
                 viennacl::matrix<T> &finalssqBeta,   //
                 viennacl::matrix<T> &finalnine,     // ssqX     
                 viennacl::matrix<T> &finalaTDa,     // ssqY    
                 viennacl::vector_base<T> &finallogD,    //of (full length)
                 viennacl::vector_base<T> &finallogP,     //(full length)    
                 viennacl::matrix<T> &finalbetahat,   //(fullsize*p, colbatch) 
                 viennacl::matrix<T> &finalLogLik,
                 const int n, 
                 const int p, 
                 const int groupsize,
                 const int colbatch,
                 const int form, // c(loglik=1, ml=2, mlFixSigma=3, mlFixBeta=4, reml=5, remlPro=6)
                 Rcpp::IntegerVector workgroupSize,
                 Rcpp::IntegerVector localSize,
                 Rcpp::IntegerVector NlocalCache,
                 const int ctx_id){
  
  
  viennacl::matrix<T> Vbatch(groupsize*n, n);
  viennacl::matrix<T> diagMat(groupsize, n);//viennacl::vector_base<T> logD(groupsize);
  viennacl::matrix<T> ab(groupsize*n, colbatch+p);             
  viennacl::matrix<T> temp2((colbatch+p)*groupsize, (colbatch+p));//viennacl::matrix<T> aTDa(groupsize, colbatch); // ssqY
  viennacl::matrix<T> temp00(groupsize*colbatch, colbatch);
  viennacl::matrix<T> temp0(groupsize, colbatch);
  viennacl::matrix<T> temp1(groupsize*n, colbatch);  // b* beta
  viennacl::matrix<T> ssqBeta0(groupsize*colbatch, colbatch); //viennacl::matrix<T> ssqBeta(groupsize, colbatch);
  viennacl::matrix<T> one(groupsize, colbatch);         
  viennacl::matrix<T> diagP(groupsize, p);
  viennacl::matrix<T> temp3(groupsize*p, colbatch);
  viennacl::matrix<T> nine0(colbatch*groupsize, colbatch); //viennacl::matrix<T> nine(groupsize, colbatch); // ssqX
  viennacl::matrix<T> two(groupsize, colbatch); //viennacl::vector_base<T> logP(groupsize);     
  viennacl::matrix<T> Qinverse(groupsize*p, p);
  viennacl::matrix<T> identity = viennacl::identity_matrix<T>(p);
  viennacl::matrix<T> QPQinverse(groupsize*p, p); //viennacl::matrix<T> betahat(groupsize*p, colbatch);  // p*groupsize   colbatch  
  viennacl::vector_base<T> logD_plusP(groupsize);     
  viennacl::matrix<T> form_temp0(groupsize, colbatch);                          
  viennacl::matrix<T> form_temp(groupsize, colbatch); 
  viennacl::matrix<T> form_temp1(groupsize, colbatch); 
  viennacl::matrix<T> form_temp2(groupsize, colbatch);
  
  Rcpp::IntegerVector localSizechol = localSize;
  localSizechol[1]= workgroupSize[1];
  
  
  const int totalsets = bigparamsBatch.size1();
  const int Niter = floor(totalsets / groupsize);
  
  viennacl::range r(0, bigparamsBatch.size2()); 
  viennacl::range r2(0, colbatch);  
  
  ///////////////////////////Loop starts !!!//////////////////////////////////////////////////////////////////////////
  for (int i=0; i< Niter; i++ ){
    
    viennacl::range ri(i*groupsize, (i+1)*groupsize);  //{0, 2}  
    viennacl::range ri_betahat(i*groupsize*p, (i+1)*groupsize*p);  //{0, 2}    
    viennacl::matrix_range<viennacl::matrix<T>> paramsBatch(bigparamsBatch, ri, r);//{m(0,0),m(0,1);m(1,0),m(1,1)} matrix_slice<matrix<double> > ms(m, s, s);         
    viennacl::matrix_range<viennacl::matrix<T>> variances(bigvariances, ri, r2);   
    viennacl::matrix_range<viennacl::matrix<T>> ssqBeta(finalssqBeta, ri, r2);    // (groupsize, colbatch)
    viennacl::matrix_range<viennacl::matrix<T>> nine(finalnine, ri, r2);  // ssqX    (groupsize, colbatch)
    viennacl::matrix_range<viennacl::matrix<T>> aTDa(finalaTDa, ri, r2);  // ssqY     (groupsize, colbatch)
    viennacl::vector_range<viennacl::vector<T> > logD(finallogD, ri);
    viennacl::vector_range<viennacl::vector<T> > logP(finallogP, ri);
    viennacl::matrix_range<viennacl::matrix<T>> betahat(finalbetahat, ri_betahat, r2);
    viennacl::matrix_range<viennacl::matrix<T>> loglik(finalLogLik, ri, r2);
    
    //42
    likfitGpu_0(Vbatch,
                coordsGpu, 
                paramsBatch, 
                diagMat, logD,
                ab,  yX, //y1,y2,y3,X               
                temp2,
                aTDa,   // ssqY
                betas, //a vclmatrix  //given by the user or provided from formula, default=null
                temp00, temp0, temp1,
                ssqBeta0, ssqBeta, //  
                one,          
                diagP,
                temp3, nine0,
                nine,      // ssqX
                two,
                logP,      
                Qinverse, identity, QPQinverse,
                betahat,  // p*rowbatch   colbatch                            
                variances,   // must be rowbatch * colbatch matrix !!!
                logD_plusP,
                form_temp0, form_temp, form_temp1, form_temp2,
                jacobian,
                loglik,
                n,  p, groupsize, colbatch,   form, // c(loglik=1, ml=2, mlFixSigma=3, mlFixBeta=4, reml=5, remlPro=6)
                workgroupSize, localSize, localSizechol, NlocalCache, ctx_id);
    
    
  }
  
  // calculate log likelihood
  
}      







//#ifdef UNDEF


//15
template<typename T> 
void likfitGpu_Templated(
    Rcpp::S4 coordsGpuR,
    Rcpp::S4 bigparamsBatchR,
    Rcpp::S4 yXR, //y1,y2,y0,X,   already transformed     
    Rcpp::S4 betasR,
    Rcpp::S4 bigvariancesR,  //already in the correct form
    Rcpp::S4 jacobianR, // already calculated, a vector 
    Rcpp::S4 finalssqBetaR,
    Rcpp::S4 finalssqXR,
    Rcpp::S4 finalssqYR,
    Rcpp::S4 finallogDR,
    Rcpp::S4 finallogPR,
    Rcpp::S4 finalbetahatR,
    Rcpp::S4 finalLogLikR,   
    const int n, 
    const int p, 
    const int groupsize,
    const int colbatch,
    const int form,
    Rcpp::IntegerVector workgroupSize,
    Rcpp::IntegerVector localSize,
    Rcpp::IntegerVector NlocalCache) {
  
  const bool BisVCL=1;
  const int ctx_id = INTEGER(finalLogLikR.slot(".context_index"))[0]-1;
  std::shared_ptr<viennacl::matrix<T> > coordsGpu = getVCLptr<T>(coordsGpuR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > bigparamsBatch = getVCLptr<T>(bigparamsBatchR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > yX = getVCLptr<T>(yXR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > betas = getVCLptr<T>(betasR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > bigvariances = getVCLptr<T>(bigvariancesR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > jacobian = getVCLptr<T>(jacobianR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > finalssqBeta = getVCLptr<T>(finalssqBetaR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > finalssqX = getVCLptr<T>(finalssqXR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > finalssqY = getVCLptr<T>(finalssqYR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::vector_base<T> > finallogD = getVCLVecptr<T>(finallogDR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::vector_base<T> > finallogP = getVCLVecptr<T>(finallogPR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > finalbetahat = getVCLptr<T>(finalbetahatR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > finalLogLik = getVCLptr<T>(finalLogLikR.slot("address"), BisVCL, ctx_id);
  
  
  likfitGpu_1<T>(*coordsGpu, 
                 *bigparamsBatch, 
                 *yX, //y1,y2,y0,X,   already transformed            
                 *betas, //a vclmatrix  //given by the user or provided from formula, default=null                          
                 *bigvariances,   // must be rowbatch * colbatch matrix !!!
                 *jacobian,
                 *finalssqBeta,
                 *finalssqX,
                 *finalssqY,
                 *finallogD,
                 *finallogP,
                 *finalbetahat,
                 *finalLogLik,
                 n, p, groupsize, colbatch, form, // c(loglik=1, ml=2, mlFixSigma=3, mlFixBeta=4, reml=5, remlPro=6)
                 workgroupSize, localSize, 
                 NlocalCache, ctx_id);
  
}










//[[Rcpp::export]]
void likfitGpu_Backend(
    Rcpp::S4 coordsGpuR,
    Rcpp::S4 bigparamsBatchR, // vclMatrix
    Rcpp::S4 yXR, //y1,y2,y0,X,   already transformed     
    Rcpp::S4 betasR,
    Rcpp::S4 bigvariancesR,  //already in the correct form
    Rcpp::S4 jacobianR, // already calculated, a vector 
    Rcpp::S4 finalssqBetaR,
    Rcpp::S4 finalssqXR,
    Rcpp::S4 finalssqYR,
    Rcpp::S4 finallogDR,
    Rcpp::S4 finallogPR,
    Rcpp::S4 finalbetahatR,
    Rcpp::S4 finalLogLikR,  
    const int n, 
    const int p, 
    const int groupsize,
    const int colbatch,
    const int form,
    Rcpp::IntegerVector workgroupSize,
    Rcpp::IntegerVector localSize,
    Rcpp::IntegerVector NlocalCache) {
  
  
  Rcpp::traits::input_parameter< std::string >::type classVarR(RCPP_GET_CLASS(finalLogLikR));
  std::string precision_type = (std::string) classVarR;
  
  
  
  if(precision_type == "fvclMatrix") {
    likfitGpu_Templated<float>(coordsGpuR, bigparamsBatchR, yXR,  betasR,  bigvariancesR,  jacobianR,    finalssqBetaR, finalssqXR,  finalssqYR,    
                               finallogDR, finallogPR, finalbetahatR, 
                               finalLogLikR, n, p, groupsize, colbatch, form, workgroupSize, localSize, 
                               NlocalCache);
  } 
  else if (precision_type == "dvclMatrix") {
    likfitGpu_Templated<double>(coordsGpuR, bigparamsBatchR, yXR,  
                                betasR,  bigvariancesR,
                                jacobianR, finalssqBetaR, finalssqXR,  finalssqYR, 
                                finallogDR, finallogPR, finalbetahatR, finalLogLikR, 
                                n, p, groupsize, colbatch, form, 
                                workgroupSize, 
                                localSize, 
                                NlocalCache);
  } else {
    Rcpp::warning("class of var must be fvclMatrix or dvclMatrix");
  }
}


//#endif
















