
#include "lgmlikFit.hpp"
//#include "gpuRandom.hpp"
//#include<Rmath.h>

using namespace Rcpp;
using namespace viennacl; 
using namespace viennacl::linalg;


//#define M_PI  3.141592653589793238463
/*
 const double PI  =3.141592653589793238463;
 const float  PI_F=3.14159265358979f;
*/
//#define M_LN_SQRT_2PI 0.918938533204672741780329736406  
// log(sqrt(2*pi)) == log(2*pi)/2


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
                 const int ctx_id, const int verbose){
  
  /* TO DO
   * compute jacobian
   */
  viennacl::vector<T> Y0 = viennacl::column(yX, 0);
  T sumLogY = viennacl::linalg::sum(viennacl::linalg::element_log(Y0));
  if(verbose) {
    Rcpp::Rcout << "\n" << "sumlog " << sumLogY << "\n" ;
    
  }
  
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
  int endThisIteration;
  
  viennacl::range r(0, bigparamsBatch.size2()); 
  viennacl::range r2(0, colbatch);  
  
  ///////////////////////////Loop starts !!!//////////////////////////////////////////////////////////////////////////
  for (int i=0; i< Niter; i++ ){
    endThisIteration = std::min((i+1)*groupsize, totalsets);   
    viennacl::range ri(i*groupsize, endThisIteration);  //{0, 2}  
    viennacl::range ri_betahat(i*groupsize*p, endThisIteration*p);  //{0, 2}    
    viennacl::matrix_range<viennacl::matrix<T>> paramsBatch(bigparamsBatch, ri, r);//{m(0,0),m(0,1);m(1,0),m(1,1)} matrix_slice<matrix<double> > ms(m, s, s);         
    viennacl::matrix_range<viennacl::matrix<T>> variances(bigvariances, ri, r2);   
    viennacl::matrix_range<viennacl::matrix<T>> ssqBeta(finalssqBeta, ri, r2);    // (groupsize, colbatch)
    viennacl::matrix_range<viennacl::matrix<T>> nine(finalnine, ri, r2);  // ssqX    (groupsize, colbatch)
    viennacl::matrix_range<viennacl::matrix<T>> aTDa(finalaTDa, ri, r2);  // ssqY     (groupsize, colbatch)
    viennacl::vector_range<viennacl::vector<T> > logD(finallogD, ri);
    viennacl::vector_range<viennacl::vector<T> > logP(finallogP, ri);
    viennacl::matrix_range<viennacl::matrix<T>> betahat(finalbetahat, ri_betahat, r2);
    viennacl::matrix_range<viennacl::matrix<T>> loglik(finalLogLik, ri, r2);
    
    
    if(verbose) {
      Rcpp::Rcout << "\n" << "iteration " << i <<
        " first parameter " << paramsBatch(0,0) << " original first parameter "
        << bigparamsBatch(0,0) << 
          "\n" ;
      
      }
    
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
    Rcpp::IntegerVector NlocalCache,
    const int verbose) {
  

  
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
                 NlocalCache, ctx_id, verbose);
  
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
    Rcpp::IntegerVector NlocalCache,
    const int verbose) {
  
  
  Rcpp::traits::input_parameter< std::string >::type classVarR(RCPP_GET_CLASS(finalLogLikR));
  std::string precision_type = (std::string) classVarR;
  
  
  
  if(precision_type == "fvclMatrix") {
    likfitGpu_Templated<float>(coordsGpuR, bigparamsBatchR, yXR,  betasR,  bigvariancesR,  jacobianR,    finalssqBetaR, finalssqXR,  finalssqYR,    
                               finallogDR, finallogPR, finalbetahatR, 
                               finalLogLikR, n, p, groupsize, colbatch, form, workgroupSize, localSize, 
                               NlocalCache, verbose);
  } 
  else if (precision_type == "dvclMatrix") {
    likfitGpu_Templated<double>(coordsGpuR, bigparamsBatchR, yXR,  
                                betasR,  bigvariancesR,
                                jacobianR, finalssqBetaR, finalssqXR,  finalssqYR, 
                                finallogDR, finallogPR, finalbetahatR, finalLogLikR, 
                                n, p, groupsize, colbatch, form, 
                                workgroupSize, 
                                localSize, 
                                NlocalCache, verbose);
  } else {
    Rcpp::warning("class of var must be fvclMatrix or dvclMatrix");
  }
}


//#endif


template <typename T> 
std::string boxcoxKernelString(int NlocalCache, int zeroCol,
                               int Nobs, int Nboxcox,
                               int Npad) {
  
  std::string typeString = openclTypeString<T>();
  std::string result = "";

  if(typeString == "double") {
    result += "\n#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n\n";
  }
  result += 
  "#define NlocalCache " + std::to_string(NlocalCache) + "\n"
  "#define DstartObs 0\n"
  "#define DstartBoxcox 0\n"
  "#define zeroCol " + std::to_string(zeroCol) + "\n"
  "#define Nobs "   + std::to_string(Nobs) + "\n"
  "#define Nboxcox "   + std::to_string(Nboxcox) + "\n"
  "#define Npad  "   + std::to_string(Npad) + "\n";
  
 result +=
   "\n__kernel void boxcox(\n"
   "__global " + typeString + " *y,\n"
   "__global " + typeString + " *param,\n"  
   "__global " + typeString + " *jacobian\n"
  "){\n\n";
   
result +=
  "int Dboxcox, Dobs;\n"
   + typeString + " logYhere, boxcoxHere;\n"
   "local " + typeString + " sumLogY[NlocalCache];\n";

 result +=
      " if(get_global_id(1)==0 & get_group_id(0) == 0){\n"
   "  sumLogY[get_local_id(0)] = 0.0;\n"
   "  for(Dobs=get_local_id(0);Dobs < Nobs; Dobs+= get_local_size(0)){\n"
   "    logYhere = log(y[DstartObs + Npad * Dobs]);\n"
   "    sumLogY[get_local_id(0)] += logYhere;\n"
   "    if(zeroCol) y[DstartObs + zeroCol + Npad * Dobs] = logYhere;\n"
   "  }// for Dobs\n"
   "  if(get_local_id(0)==0){\n"
   "    for(Dobs=1;Dobs < get_local_size(0);Dobs++){\n"
   "     sumLogY[0] += sumLogY[Dobs];\n"
   "    }// for Dobs\n"
   "    if(get_global_id(0)==0){\n"
   "     for(Dboxcox=0;Dboxcox < Nboxcox;Dboxcox++){\n"
   " //   jacobian = -2*(BoxCox-1)* sum(log(yXcpu[,1]))\n"
   "      jacobian[DstartBoxcox + Dboxcox] =  -2*(param[DstartBoxcox + Dboxcox]-1)*sumLogY[0];\n"
   "     }// for Dboxcox\n"
   "    }// if global0\n"
   "  }// if local0\n"
   "}// if global and group\n";
   
   result +=
   "for(Dboxcox = get_global_id(1)+1; Dboxcox < Nboxcox; Dboxcox+= get_global_size(1)){\n"
   " boxcoxHere = param[DstartBoxcox + Dboxcox];\n"
   " if(Dboxcox != zeroCol) {\n"
   "  for(Dobs=get_global_id(0);Dobs < Nobs; Dobs+= get_global_size(0)){\n"
   "//    transformed_y[ ,i] <- ((yXcpu[ ,1]^BoxCox[i]) - 1)/BoxCox[i]\n"
   "      y[DstartObs + Dboxcox + Npad * Dobs] = \n"
   "        (pow(y[DstartObs + Npad * Dobs],boxcoxHere)-1)/boxcoxHere;\n"
   "    }// for Dobs\n"
   "  }// if zerocol\n"
   "}// for Dboxcox\n";

  result +=   "}// kernel\n";
  return(result);
}


/*
 * sum logs of rows
 * use only one work item dimension
 * groups are rows, local items sum over colums
 * 
 */
template <typename T> 
std::string logRowSumString(int NlocalCache) {  
  
  std::string typeString = openclTypeString<T>();
  
  
  std::string result = "";
  
  if(typeString == "double") {
    result += "\n#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n";
  }
  result += 
    "#define NlocalCache " + std::to_string(NlocalCache) + "\n";    

  result += 
    "\n\n__kernel void logRowSum(\n"
    "  __global " + typeString + "* x,\n"  
    "  __global " + typeString + "* output,\n" 
    "  int Nrow, int Ncol, int Npad, int startX, int startOutput"
    "){\n\n";  
  
  result += "int Drow, Dcol, Dindex;\n";
  result += "local " + typeString + " thesum[NlocalCache];\n";
  
  result += 
    
    "  for(Drow = get_group_id(0);   Drow < Nrow;    Drow+=get_num_groups(0)){\n"
    "    thesum[get_local_id(0)] = 0.0;\n"
    "    for(Dcol = get_global_id(0),   Dindex = startX + Drow*NpadCol+Dcol;\n" 
    "        Dcol < Ncol; Dcol+=get_global_size(0), Dindex+=get_global_size(0)){\n"
    "        thesum[get_local_id(0)] += log(x[Dindex]);\n"
    "    } // end loop through columns\n"
    "  if(get_local_id(0) == 0) {\n"
    "    for(Dcol = 1; Dcol < get_local_size(0); Dcol++){\n"
    "      thesum[0] += thesum[Dcol];\n"
    "    }//for\n"
    "    output[Drow + startOutput] = thesum[0];\n"
    "  }//if\n"
    "  } // end loop through rows\n"
    "}\n";
  
  return(result);
}

template<typename T> 
void addBoxcoxToData(
    viennacl::matrix_base<T> &yx,
    viennacl::vector_base<T>  &boxcox,
    viennacl::vector_base<T> &jacobian,
    Rcpp::IntegerVector workgroupSize, 
    Rcpp::IntegerVector localSize, 
    Rcpp::IntegerVector NlocalCache, 
    const int ctx_id,
    Rcpp::IntegerVector verbose){

  viennacl::ocl::switch_context(ctx_id);
  viennacl::ocl::context ctx(viennacl::ocl::get_context(ctx_id));
  const int Ndatasets = boxcox.size();
  if(verbose[0] & boxcox.size() > 1) {
    if(boxcox[1] != 0){
      Rcpp::warning("second entry of boxcox parameters should be zero");
    }
  }
  
  std::string theBoxcoxKernel = boxcoxKernelString<T>(
    NlocalCache[0],
               1, yx.size1(), Ndatasets, 
               yx.internal_size2());
  
  if(verbose[0]>1) {
    Rcpp::Rcout << "\n" << theBoxcoxKernel << "\n";
  }
  viennacl::ocl::program & my_prog_boxcox = viennacl::ocl::current_context().add_program(theBoxcoxKernel, "mkb");
  
  viennacl::ocl::kernel & boxcoxKernel = my_prog_boxcox.get_kernel("boxcox");
  boxcoxKernel.global_work_size(0, (cl_uint) (workgroupSize[0] ) );
  boxcoxKernel.global_work_size(1, (cl_uint) (workgroupSize[1] ) );
  boxcoxKernel.local_work_size(0, (cl_uint) (localSize[0]));
  boxcoxKernel.local_work_size(1, (cl_uint) (localSize[1]));

  viennacl::ocl::enqueue(
    boxcoxKernel(
      yx, boxcox, jacobian));
  
  }


template<typename T> 
void likfitGpuP(viennacl::matrix_base<T> &yx, 
                viennacl::matrix_base<T> &coords, 
                viennacl::matrix_base<T> &params, 
                  viennacl::matrix_base<T> &betas,
                  viennacl::matrix_base<T> &ssqY,
                  viennacl::matrix_base<T> &ssqX,
                  viennacl::vector_base<T> &detVar,
                  viennacl::matrix_base<T> &detReml,
                  int Ndatasets,
                  Rcpp::IntegerVector NparamPerIter,
                  Rcpp::IntegerVector workgroupSize, 
                  Rcpp::IntegerVector localSize, 
                  Rcpp::IntegerVector NlocalCache, 
                  const int ctx_id, 
                  Rcpp::IntegerVector verbose){
     int Nobs = yx.size1();
   int Nparams = params.size1();
   int Ncovariates = yx.size2() - Ndatasets;

     int Niter = ceil(Nparams / NparamPerIter[0]);
    int Diter;
    int endThisIteration;
    int DiterIndex, NthisIteration;
    int verboseMatern = verbose[0]>1;
    
    viennacl::ocl::switch_context(ctx_id);
    viennacl::ocl::context ctx(viennacl::ocl::get_context(ctx_id));

    viennacl::matrix<T> Vbatch(NparamPerIter[0]*Nobs, Nobs);
    viennacl::matrix<T> cholDiagMat(NparamPerIter[0], Nobs);
    viennacl::matrix<T> LinvYX(NparamPerIter[0]*Nobs, yx.size2());
    viennacl::matrix<T> ssqYX(NparamPerIter[0]*yx.size2(), yx.size2());
    
    fill22params(params);
    
    
  /* 
   * compute boxcox and jacobian
   */
  
  // create and compile kernels
   int Ncell = Nobs * (Nobs - 1)/2, maxIter = 1500;

  std::string maternClString = maternBatchKernelString<T>(
    maxIter,
    Nobs, Ncell, 
    NparamPerIter[0],//NmatrixMax
    Vbatch.internal_size2(), //NpadVbatch
    Vbatch.internal_size2()*Nobs, //NpadBetweenMatrices,
    coords.internal_size2(), 
    params.internal_size2(), //NpadParams,
    localSize[0],
    NlocalCache[0],
    1L, 1L, 1L, 0L
  );
  int allowOverflow = ( ((int) Vbatch.size2() ) > NlocalCache[0] );
  std::string cholClString = cholBatchKernelString<T>(
    0L, // colstart
    Nobs, // colend
    Nobs, // N
    Vbatch.internal_size2(), // Npad
    cholDiagMat.internal_size2(), // NpadDiag
    Vbatch.size2() * Vbatch.internal_size2(),// NpadBetweenMatrices,
    0L,//NstartA
    0L,// NstartD
    NlocalCache, //Ncache
    localSize, //Nlocal
    allowOverflow, // allowoverflow
    1L // do log determinant
  );

  const int Ncol = yx.size2();
  const int Ngroups1 = static_cast<T>(workgroupSize[1]) / static_cast<T>(localSize[1]);
  const int NcolsPerGroup = std::ceil( static_cast<T>(Ncol) / static_cast<T>(Ngroups1));
  const int NrowsToCache = std::floor(static_cast<T>(NlocalCache[0]) /static_cast<T>(NcolsPerGroup));
  
  std::string backsolveString = backsolveBatchString<T>(
      0,// sameB,
      1,// diagIsOne,
      Nobs,// Nrow, 
      Ncol,
      LinvYX.internal_size2(),// NpadC, 
      Vbatch.internal_size2(),// NpadA, 
      yx.internal_size2(),// NpadB,
      LinvYX.internal_size2()*Nobs,// NpadBetweenMatricesC,
      Vbatch.internal_size2()*Nobs,// NpadBetweenMatricesA,
      yx.internal_size2()*Nobs,// NpadBetweenMatricesB,
      0,// NstartC,
      0,// NstartA,
      0,// NstartB,
      NrowsToCache, 
      NcolsPerGroup,
      NcolsPerGroup * NrowsToCache,// NlocalCacheC,  
      localSize[0] * localSize[1] * NcolsPerGroup,// NlocalCacheSum,   
      localSize[0] * localSize[1]//NpadBetweenMatricesSum 
  );
    
std::string crossprodKernelString = crossprodBatchString<T>(
  Nobs,//const int Nrow, 
  yx.size2(),//const int Ncol,
  //    const int Nmatrix, 
  ssqYX.internal_size2(),//const int NpadC, 
  LinvYX.internal_size2(),//const int NpadA,
  cholDiagMat.internal_size2(),//const int NpadD, // set to zero to omit D
  1, //const int invertD, // set to 1 for A^T D^(-1) A
  0,//const int NstartC,  // newly added
  0,//const int NstartA,  // new
  0,//const int NstartD,  // new
  ssqYX.internal_size2()*Nobs, //const int NpadBetweenMatricesC,
  LinvYX.internal_size2()*Nobs, //const int NpadBetweenMatricesA,
  NlocalCache[0]/2, // NlocalCacheA, 
  localSize// Nlocal// cache a Nlocal[0] by Nlocal[1] submatrix of C
);
        
  if(verbose[0]>1) {
      Rcpp::Rcout << maternClString << "\n";
    Rcpp::Rcout << cholClString << "\n";
  }
  
  
  viennacl::ocl::program & my_prog_matern = viennacl::ocl::current_context().add_program(maternClString, "mykernelmatern");
  viennacl::ocl::program & my_prog_chol = viennacl::ocl::current_context().add_program(cholClString, "mykernelchol");
  viennacl::ocl::program & my_prog_backsolve = viennacl::ocl::current_context().add_program(backsolveString, "mykernelbacksolve");
  viennacl::ocl::program & my_prog_crossprod = viennacl::ocl::current_context().add_program(crossprodKernelString, "mykernelcrossprod");
    
  viennacl::ocl::kernel & maternKernel = my_prog_matern.get_kernel("maternBatch");
  viennacl::ocl::kernel & cholKernel = my_prog_chol.get_kernel("cholBatch");
  viennacl::ocl::kernel & backsolveKernel = my_prog_backsolve.get_kernel("backsolveBatch");
  viennacl::ocl::kernel & crossprodKernel = my_prog_crossprod.get_kernel("crossprodBatch");
  
  
  // dimension 0 is cell, dimension 1 is matrix
  maternKernel.global_work_size(0, workgroupSize[0] ); 
  maternKernel.global_work_size(1, workgroupSize[1] ); 
  maternKernel.local_work_size(0, localSize[0]);
  maternKernel.local_work_size(1, localSize[1]);
  cholKernel.global_work_size(0, workgroupSize[0] ); 
  cholKernel.global_work_size(1, workgroupSize[1] ); 
  cholKernel.local_work_size(0, localSize[0]);
  cholKernel.local_work_size(1, localSize[1]);
  backsolveKernel.global_work_size(0, workgroupSize[0] ); 
  backsolveKernel.global_work_size(1, workgroupSize[1] ); 
  backsolveKernel.local_work_size(0, localSize[0]);
  backsolveKernel.local_work_size(1, localSize[1]);
  crossprodKernel.global_work_size(0, workgroupSize[0] ); 
  crossprodKernel.global_work_size(1, workgroupSize[1] ); 
  crossprodKernel.local_work_size(0, localSize[0]);
  crossprodKernel.local_work_size(1, localSize[1]);
  
  
  ///////////////////////////Loop starts !!!//////////////////////////////////////////////////////////////////////////
  for (Diter=0,DiterIndex=0; Diter< Niter; Diter++,DiterIndex += NparamPerIter[0]){

    endThisIteration = std::min(DiterIndex + NparamPerIter[0], Nparams);
    NthisIteration = endThisIteration - DiterIndex;

    if(verbose[0]) {
      Rcpp::Rcout << "\n" << "Diter " << Diter <<" DiterIndex " << DiterIndex << " endThisIteration " << 
        endThisIteration << " Nthisiteration " << NthisIteration <<"\n";
    }
    
    
    viennacl::ocl::enqueue(maternKernel(Vbatch, coords, params, DiterIndex, NthisIteration));

    //#Vbatch=LDL^T, cholesky decomposition
    viennacl::ocl::enqueue(cholKernel(Vbatch, cholDiagMat, NthisIteration, detVar, DiterIndex));
    // LinvYX = L^(-1) YX
    viennacl::ocl::enqueue(backsolveKernel(LinvYX, Vbatch, yx, NthisIteration));
    // ssqYX = YX^Y L^(-1)T D^(-1) L^(-1) YX
    viennacl::ocl::enqueue(crossprodKernel(ssqYX, LinvYX, cholDiagMat, NthisIteration));
    
    // copy diag(ssqYX[1:Ndatasets, 1:Ndatasets]) to ssqY
    
    
    } // Diter
}

template<typename T> 
void likfitGpuP_Templated(
    Rcpp::S4 yx,
    Rcpp::S4 coords,
    Rcpp::S4 params,
    Rcpp::S4 boxcox,
    Rcpp::S4 betas,
    Rcpp::S4 ssqY,
    Rcpp::S4 ssqX,
    Rcpp::S4 detVar,
    Rcpp::S4 detReml,
    Rcpp::S4 jacobian,
    Rcpp::IntegerVector NparamPerIter,
    Rcpp::IntegerVector workgroupSize,
    Rcpp::IntegerVector localSize,
    Rcpp::IntegerVector NlocalCache,
    Rcpp::IntegerVector verbose) {
  
  
  const bool BisVCL=1;
  const int ctx_id = INTEGER(yx.slot(".context_index"))[0]-1;
  std::shared_ptr<viennacl::matrix<T> > yxGpu = getVCLptr<T>(yx.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > coordsGpu = getVCLptr<T>(coords.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > paramsGpu = getVCLptr<T>(params.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > betasGpu = getVCLptr<T>(betas.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > ssqYGpu = getVCLptr<T>(ssqY.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > ssqXGpu = getVCLptr<T>(ssqX.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::vector_base<T> > detVarGpu = getVCLVecptr<T>(detVar.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > detRemlGpu = getVCLptr<T>(detReml.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::vector_base<T> > boxcoxGpu = getVCLVecptr<T>(boxcox.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::vector_base<T> > jacobianGpu = getVCLVecptr<T>(jacobian.slot("address"), BisVCL, ctx_id);

  addBoxcoxToData<T>(
    *yxGpu,
    *boxcoxGpu,
    *jacobianGpu,
    workgroupSize, 
    localSize,
    NlocalCache,
    ctx_id,
    verbose);
    
  likfitGpuP<T>(
                 *yxGpu, 
                 *coordsGpu, 
                 *paramsGpu, 
                 *betasGpu,
                 *ssqYGpu, *ssqXGpu,
                 *detVarGpu, *detRemlGpu,
                 (*boxcoxGpu).size(),// Ndatasets
                 NparamPerIter,
                 workgroupSize, 
                 localSize, 
                 NlocalCache, ctx_id, verbose);
  
  }






//[[Rcpp::export]]
void likfitGpu_BackendP(
    Rcpp::S4 yx,
    Rcpp::S4 coords,
    Rcpp::S4 params,
    Rcpp::S4 boxcox,
    Rcpp::S4 betas,
    Rcpp::S4 ssqY,
    Rcpp::S4 ssqX,
    Rcpp::S4 detVar,
    Rcpp::S4 detReml,
    Rcpp::S4 jacobian,
    Rcpp::IntegerVector NparamPerIter,
    Rcpp::IntegerVector workgroupSize,
    Rcpp::IntegerVector localSize,
    Rcpp::IntegerVector NlocalCache,
    Rcpp::IntegerVector verbose) {
  
  
  Rcpp::traits::input_parameter< std::string >::type classVarR(RCPP_GET_CLASS(yx));
  std::string precision_type = (std::string) classVarR;
  
  
  if(precision_type == "fvclMatrix") {
    likfitGpuP_Templated<float>(  yx,
                                  coords,
                                  params,
                                  boxcox,
                                  betas,
                                  ssqY,
                                  ssqX,
                                  detVar,
                                  detReml,
                                  jacobian,
                                  NparamPerIter,
                                  workgroupSize,
                                  localSize,
                                  NlocalCache,
                                  verbose);
  } 
  else if (precision_type == "dvclMatrix") {
    likfitGpuP_Templated<double>(  yx,
                                   coords,
                                   params,
                                   boxcox,
                                   betas,
                                   ssqY,
                                   ssqX,
                                   detVar,
                                   detReml,
                                   jacobian,
                                   NparamPerIter,
                                   workgroupSize,
                                   localSize,
                                   NlocalCache,
                                   verbose);
  } else {
    Rcpp::warning("class of var must be fvclMatrix or dvclMatrix");
  }
}


//#endif
















