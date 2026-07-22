/* type = 0, distances is a vector
 type = 1, distances is a symmetric matrix,
 only lower triangle is used
 type = 2, as type=1 but return the cholesky
 type = 3, as type=1, but return inverse
 nugget is only used if type >1
 and added to the diagonals

If the cholkesy is performed (type > 1):
 halfLogDet is the log determinant if the cholesky matrix
 type is info from dpotrf
if the precision is computed type is info from dpotrfi
... and if chol of precision is computed, type is from dtrtri
 */

#include"geostatsp.h"

void maternArasterBpoints(
    double *Axmin, double *Axres,
    int *AxN,
    double *Aymax, double *Ayres, int *AyN,
    double *Bx, double *By, int *BN,
    double *result,
    double  *range, double*shape, double *variance,
    double *anisoRatio, double *anisoAngleRadians) {

  int DB, DAx, DAy, AyN2, AxN2, BN2;
  int Dindex,Ncell;
  double distCellRight[2], distCellDown[2], distTopLeft[2], distRowHead[2];
  double distTopLeftR[2], distHere[2];
  double costheta, sintheta, anisoRatioSq, distSq;
  maternParams mp;

  maternParamsSet(&mp, *range, *shape, *variance);

  AyN2 = *AyN;
  AxN2 = *AxN;
  Ncell = AyN2*AxN2;
  BN2 = *BN;
  *Axmin += *Axres/2; // add half a cells size, assign each cell
  *Aymax -= *Ayres/2; // it's centroid rather than it's top left corner


  costheta = cos(*anisoAngleRadians);
  sintheta = sin(*anisoAngleRadians);
  anisoRatioSq = (*anisoRatio)*(*anisoRatio);

  distCellRight[0] = costheta *(*Axres);
  distCellRight[1] = sintheta * (*Axres);

  distCellDown[0] =  sintheta * (*Ayres);
  distCellDown[1] =  - costheta * (*Ayres);

  for(DB=0;DB<BN2;++DB){ // loop through points
      Dindex = DB*Ncell;
      distTopLeft[0]= (Bx[DB]-*Axmin); // distance from point DB to
      distTopLeft[1]= (By[DB]-*Aymax); //   top left corner of raster

      distTopLeftR[0]= costheta * distTopLeft[0] - sintheta * distTopLeft[1];
      distTopLeftR[1] = sintheta *distTopLeft[0] + costheta * distTopLeft[1];

      distRowHead[0] = distTopLeftR[0]; // distance to leftmost cell of row DAy
      distRowHead[1] = distTopLeftR[1];

      for(DAy=0;DAy<AyN2;++DAy){ // loop through y of raster

          distHere[0] = distRowHead[0]; // dist to cell DAx Day
          distHere[1] = distRowHead[1];
          for(DAx=0;DAx<AxN2;++DAx){ // loop through x of raster

              distSq = distHere[0]*distHere[0] +
                  distHere[1]*distHere[1]/anisoRatioSq;
              result[Dindex] = maternPoint(distSq, &mp);

              ++Dindex;
              distHere[0] -= distCellRight[0];
              distHere[1] -= distCellRight[1];
          }

          distRowHead[0] -=distCellDown[0];
          distRowHead[1] -=distCellDown[1];

      }

  }

}

// returns an N by N matrix for matern correlation
// for N points with vectors of coordinates x and y
// type = 0 or 1 return correlation,
// type=2 chol, type=3 precision, type 4 chol of precsion
void maternAniso(
    const double *x,
    const double *y,
    const int *N,
    double *result,
    const double *range,
    const double *shape,
    const double *variance,
    const double *anisoRatio,
    const double *anisoAngleRadians,
    const double *nugget,
    int *type,
    double *halfLogDet
) {
  // type=2 return cholesky
  int Drow, Dcol, Nm1, Dcolp1, N2;
  int Dindex;

  double anisoRatioSq, dist[2], distRotate[2], costheta, sintheta, distSq;
  maternParams mp;

  maternParamsSet(&mp, *range, *shape, *variance);

  costheta = cos(*anisoAngleRadians);
  sintheta = sin(*anisoAngleRadians);
  anisoRatioSq = (*anisoRatio)*(*anisoRatio);

  Nm1 = *N-1;
  N2 = *N;

  result[N2*N2-1] = *variance + *nugget;
  for(Dcol=0;Dcol < Nm1;++Dcol) {
      Dcolp1 = Dcol + 1;
      Dindex = Dcol*N2+Dcol;
      result[Dindex] = *variance + *nugget;
      for(Drow=Dcolp1; Drow < N2; ++Drow) {
          Dindex++;
          dist[0] = x[Dcol] - x[Drow];
          dist[1] = y[Dcol] - y[Drow];

          // rotate anticlockwise by aniso.angle.radians
          // distRotate =  ( cos(theta)  -sin(theta)  )  dist
          //               ( sin(theta)   cos(theta)  )
          distRotate[0] = costheta *dist[0] - sintheta * dist[1];
          distRotate[1] = sintheta *dist[0] + costheta * dist[1];

          distSq = distRotate[0]*distRotate[0] +
              distRotate[1]*distRotate[1]/anisoRatioSq;
          result[Dindex] = maternPoint(distSq, &mp);
      } // end for Drow
  }// end for Dcol

  if(*type >1 ){ // cholesky
      int requestedType = *type;
      int infoLapack;
      F77_CALL(dpotrf)("L", N, result, N, &infoLapack FCONE);
      if(infoLapack != 0){
          *type = infoLapack;
          *halfLogDet = NA_REAL;
          return;
      }
      *halfLogDet=0;  // the log determinant

      for(Drow = 0; Drow < N2; Drow++)
        *halfLogDet += log(result[Drow*N2+Drow]);
      if(requestedType == 3){ // precision
          F77_NAME(dpotri)("L", N,
              result, N, &infoLapack FCONE);
          if(infoLapack != 0){
              *type = infoLapack;
              *halfLogDet = NA_REAL;
              return;
          }
      } else if (requestedType==4) {// cholkesy of precision
            F77_NAME(dtrtri)("L", "N",N,
              result, N, &infoLapack FCONE FCONE);
            if(infoLapack != 0){
                *type = infoLapack;
                *halfLogDet = NA_REAL;
                return;
            }
      } else {
          Drow = 0;
      }
      *type = 0;
  }

//  free(bk);
}

void matern(
    const double *distance,
    const int *N,
    double *result,
    const double *range,
    const double *shape,
    const double *variance,
    const double *nugget,
    int *type,
    double *halfLogDet) {
  int D, Dcol, Ncol, Nrow, rowEnd, addToRowStart;
  double distAbs;
  maternParams mp;

  /* evaluate the matern:
     thisx = abs(x)*(sqrt(8*shape)/range)
     result = variance/(gamma(shape)*2^(shape-1)) * thisx^shape * K(thisx, shape)
   */
  maternParamsSet(&mp, *range, *shape, *variance);

  Nrow = *N;
  if(*type){ // lower triangle
      for(Dcol=0;Dcol<Nrow;++Dcol){
          // diagonals
          result[Dcol*Nrow+Dcol] =
              *variance + *nugget;
      }
      Ncol = *N-1; // the last column isn't done
      addToRowStart = 1; // the diagonals aren't done with materns
  } else {
      Ncol=1;
      addToRowStart=0;
  }

  for(Dcol=0;Dcol<Ncol;++Dcol) {
      rowEnd = Nrow*Dcol+Nrow;
      for(D=Dcol*Nrow+Dcol+addToRowStart; D < rowEnd; D++) {
          distAbs = fabs(distance[D]);
          result[D] = maternPoint(distAbs * distAbs, &mp);
      } //D
  } // Dcol

  if(*type >1 ){ // cholesky
      int requestedType = *type;
      int infoLapack;
      F77_CALL(dpotrf)("L", N, result, N, &infoLapack FCONE);
      if(infoLapack != 0){
          *type = infoLapack;
          *halfLogDet = NA_REAL;
          return;
      }
      *halfLogDet=0;  // the log determinant
      for(D = 0; D < Nrow; D++)
        *halfLogDet += log(result[D*Nrow+D]);
      if(requestedType == 3){//precision
          F77_NAME(dpotri)("L", N,
              result, N, &infoLapack FCONE);
          if(infoLapack != 0){
              *type = infoLapack;
              *halfLogDet = NA_REAL;
              return;
          }
      } else if (requestedType==4) {// cholesky of precision
          F77_NAME(dtrtri)("L", "N",N,
              result, N, &infoLapack FCONE FCONE);
          if(infoLapack != 0){
              *type = infoLapack;
              *halfLogDet = NA_REAL;
              return;
          }
      } else {
          D = 0;
      }
      *type = 0;
  }
  //free(bk);
}

// matern for a vector of parameters
void maternForL(
    const double *xcoord,
    const double *ycoord,
    const int *N,
    double *corMat,
    const double *param,
    // nugget, variance,
    // range, shape,
    // anisoRatio, ansioAngleRadians
    const int *aniso,
    const int *withoutNugget,
    int *type,
    double *halfLogDet
){


  double nugget;

  if(*withoutNugget){
      nugget = 0.0;
  } else {
      nugget = param[0];
  }

  if(*aniso) {
      maternAniso(
          xcoord,ycoord,
          N,
          corMat,
          &param[2],
          &param[3],
          &param[1],
          &param[4],
          &param[5],
          &nugget,
          type,
          halfLogDet);
  } else {
      matern(xcoord,
             N,
             corMat,
             &param[2],
             &param[3],
             &param[1],
             &nugget,
             type,
             halfLogDet);
  }


}

int typeStringToInt(SEXP type){

  const char *typeSeq[] = {"variance", "cholesky", "precision","inverseCholesky"};
  int D, typeInt = -1;

  // find type integer, and check it's valid
  for(D = 0;D< 4;++D){
      if( strcmp(CHAR(STRING_ELT(type, 0)), typeSeq[D]) == 0 ) {
          typeInt = D;
      }
  }
  if(typeInt>3) {
      Rprintf("d %s %s", CHAR(STRING_ELT(type, 0)), typeSeq[0]);
      error("'type' in matern, must be one of 'variance','cholesky','precision','inverseCholesky'");
  } else {
      // indexed from 1
      typeInt = typeInt + 1;
  }
  return typeInt;
}



SEXP maternPoints(
    SEXP points,
    SEXP result,
    SEXP param,
    // range,
    // shape,
    // variance,
    // anisoRatio,
    // anisoAngleRadians,
    // nugget,
    SEXP type) {

  SEXP halfLogDet= PROTECT(NEW_NUMERIC(1));
  int N= Rf_nrows(points);
  int requestedType = INTEGER(type)[0];

  
  

  maternAniso(
      REAL(points), // x
      &REAL(points)[N], // y
      &N,
      REAL(GET_SLOT(result, install("x"))),
      &REAL(param)[0],// range,
      &REAL(param)[1],// shape,
      &REAL(param)[2],// variance,
      &REAL(param)[3],// anisoRatio,
      &REAL(param)[4],// anisoAngleRadians,
      &REAL(param)[5],// nugget,
      INTEGER(type),
      REAL(halfLogDet)
  );
  if(requestedType > 1 && INTEGER(type)[0] != 0){
      error("LAPACK factorization failed in maternPoints, info=%d", INTEGER(type)[0]);
  }
  UNPROTECT(1);
  return halfLogDet;
}


SEXP maternDistance(
    SEXP distance,
    SEXP result,
    SEXP param,
    // range, shape,
    // variance, nugget,
    int *type
    //c('variance','cholesky','precision','inverseCholesky')
) {

//  SEXP halfLogDet= PROTECT(NEW_NUMERIC(1));
  double halfLogDet;
  const char
  *valid[] = {"dsyMatrix"};
  int typeInt=*type, requestedType=*type, D, D2, N;
  double *P, *Presult;


  
  N = INTEGER(GET_SLOT(distance, install("Dim")))[0];
  P = REAL(GET_SLOT(distance, install("x")));
  Presult = REAL(GET_SLOT(result, install("x")));

  // check distance is symmetric
  if(R_check_class_etc(distance, valid)) {
      error("invalid class of 'distance' in maternDistance, must be dsyMatrix");
  }

  // check values stored in lower triangle
  if(strcmp(
      CHAR(STRING_ELT(GET_SLOT(distance, install("uplo")), 0) ),
      "L"
  ) != 0) {
      // not lower triangle, copy over
      for (D = 1; D < N; D++) // rows
        for (D2 = 0; D2 < D && D2 < N; D2++) // columns
          P[D + D2*N] = P[D2 + D*N];
  }




  matern(
      P,
      &N, //N
      Presult,
      &REAL(param)[0],// range,
      &REAL(param)[1],// shape,
      &REAL(param)[2],// variance,
      &REAL(param)[3],// nugget,
      &typeInt,
      &halfLogDet //REAL(halfLogDet)
      );
  if(requestedType > 1 && typeInt != 0){
      error("LAPACK factorization failed in maternDistance, info=%d", typeInt);
  }
//  UNPROTECT(1);
  return Rf_ScalarReal(halfLogDet);
}


void maternRaster(
    double *Axmin, 
    double *Axres,
    int *AxN,
    double *Aymax, 
    double *Ayres, 
    int *AyN,
    double *result,
    double  *range, 
    double *shape, 
    double *variance,
    double *anisoRatio, 
    double *anisoAngleRadians,
    int *type) {

  int DB, DBx, DBy, DAx, DAy, AyN2, AxN2;
  int Dindex,Ncell;
  double distCellRight[2], distCellDown[2], distTopLeft[2], distRowHead[2];
  double distTopLeftR[2], distHere[2], Bx, By;
  double costheta, sintheta, anisoRatioSq, distSq;
  maternParams mp;

  maternParamsSet(&mp, *range, *shape, *variance);

  AyN2 = *AyN;
  AxN2 = *AxN;
  Ncell = AyN2*AxN2;
  *Axmin += *Axres/2; // add half a cells size, assign each cell
  *Aymax -= *Ayres/2; // it's centroid rather than it's top left corner


  costheta = cos(*anisoAngleRadians);
  sintheta = sin(*anisoAngleRadians);
  anisoRatioSq = (*anisoRatio)*(*anisoRatio);

  distCellRight[0] = costheta *(*Axres);
  distCellRight[1] = sintheta * (*Axres);

  distCellDown[0] =  sintheta * (*Ayres);
  distCellDown[1] =  - costheta * (*Ayres);

  DB = 0; // cell index
  for(DBy=0;DBy<AyN2;++DBy){ // loop through rows of raster B
      By = *Aymax - DBy * (*Ayres);
      for(DBx=0;DBx<AxN2;++DBx){ // loop through columns of raster B
          Bx = *Axmin + DBx * (*Axres);
          distTopLeft[0]= (Bx-*Axmin); // distance from point DB to
          distTopLeft[1]= (By-*Aymax); //   top left corner of raster

          distTopLeftR[0]= costheta * distTopLeft[0] - sintheta * distTopLeft[1];
          distTopLeftR[1] = sintheta *distTopLeft[0] + costheta * distTopLeft[1];

          // distance to leftmost cell of row DAy

          distRowHead[0] = distTopLeftR[0] - distCellDown[0]*DBy;
          distRowHead[1] = distTopLeftR[1] - distCellDown[1]*DBy;

          for(DAy=DBy;DAy<AyN2;++DAy){ // loop through y of raster A
              Dindex = DB*Ncell + DAy*AxN2; // index in covariance matrix

              distHere[0] = distRowHead[0]; // dist to cell DAx Day
              distHere[1] = distRowHead[1];
              for(DAx=0;DAx<AxN2;++DAx){ // loop through x of raster A

                  distSq = distHere[0]*distHere[0] +
                      distHere[1]*distHere[1]/anisoRatioSq;
                  result[Dindex] = maternPoint(distSq, &mp);

                  ++Dindex;
                  distHere[0] -= distCellRight[0];
                  distHere[1] -= distCellRight[1];
              } // x of raster A

              distRowHead[0] -=distCellDown[0];
              distRowHead[1] -=distCellDown[1];

          } // loop through y of raster A
          DB += 1;
      } // loop through x of raster B
  } // y of raster B

  if(*type >1 ){ // cholesky
      int requestedType = *type;
      int infoLapack;
      F77_CALL(dpotrf)("L", &Ncell, result, &Ncell, &infoLapack FCONE);
      if(infoLapack != 0){
          *type = infoLapack;
          error("LAPACK dpotrf failed in maternRaster, info=%d", infoLapack);
      }
      if(requestedType == 3){//precision
          F77_NAME(dpotri)("L", &Ncell,
              result, &Ncell,&infoLapack FCONE);
          if(infoLapack != 0){
              *type = infoLapack;
              error("LAPACK dpotri failed in maternRaster, info=%d", infoLapack);
          }
      } else if (requestedType==4) {// cholesky of precision
          F77_NAME(dtrtri)("L", "N", &Ncell,
              result, &Ncell,&infoLapack FCONE FCONE);
          if(infoLapack != 0){
              *type = infoLapack;
              error("LAPACK dtrtri failed in maternRaster, info=%d", infoLapack);
          }
      }
      *type = 0;
  }

}


void maternRasterConditional(
    double *Axmin, double *Axres, int *AxN,
    double *Aymax, double *Ayres, int *AyN,
    double *ydata, double *yx, double *yy, int *Ny,
    double *result,
    int *Nsim, // number of realisations per parameter set
    int *Nparam, // number of parameter sets
    double *nugget,
    double  *range, double*shape, double *variance,
    double *anisoRatio, double *anisoAngleRadians,
    double *inVarGrid
) {

  int oneI=1, fourI=4, Ncell, Nrandom;
  int Dparam, D;
  double *resultHere, *ydataHere, oneD=1.0, minusOneD=-1.0;
  //double *varY, *covDataGrid
  double *varGrid, halfLogDet=0.0;

  Ncell = (*AyN) * (*AxN);
  Nrandom = Ncell * (*Nsim);
//  varY = (double *) calloc((*Ny)*(*Ny), sizeof(double));
//  covDataGrid = (double *) calloc( (*Ny) * Ncell, sizeof(double));
  SEXP varY = PROTECT(NEW_NUMERIC( (*Ny)*(*Ny) ));
  SEXP covDataGrid = PROTECT(NEW_NUMERIC((*Ny) * Ncell));

  varGrid = inVarGrid;//(double *) calloc(NcellSq, sizeof(double));


  for(Dparam=0; Dparam < *Nparam; ++Dparam) {
      resultHere = &result[Dparam*Nrandom];
      ydataHere = &ydata[Dparam* (*Ny)];

      // random noise
      for(D=0;D<Nrandom;++D) {
          resultHere[D] = norm_rand();
      }
      // var Y

      maternAniso(
          yx, yy, Ny,
          REAL(varY),
          &range[Dparam],
          &shape[Dparam],
          &variance[Dparam],
          &anisoRatio[Dparam],
          &anisoAngleRadians[Dparam],
          &nugget[Dparam], &fourI,
          &halfLogDet
      );

      // covUY
      maternArasterBpoints(
          Axmin, Axres, AxN, Aymax, Ayres, AyN,
          yx, yy, Ny,
          REAL(covDataGrid),
          &range[Dparam], &shape[Dparam], &variance[Dparam],
          &anisoRatio[Dparam], &anisoAngleRadians[Dparam]);

      /*
       * varY = L Lt, varY^(-1) = Lt^(-1) L^(-1)
       * want Linv %*% t(covUV) or t( covUV %*% t(Linv) )

       * side = 'R' uplo='L' transa='T' diag='N'
       * M = Ngrid N=Ny
       * alpha = 1.0
       * A= varY (really Linv)
       * LDA=Ny
       * B = covDataGrid
       * LDB = Ngrid

    DTRMM  performs one of the matrix-matrix operations

        B := alpha*op( A )*B,   or   B := alpha*B*op( A ),

       */

      F77_NAME(dtrmm)(
          "R", "L", "T", "N",
          &Ncell, Ny, &oneD,
          REAL(varY), Ny,
          REAL(covDataGrid), &Ncell
          FCONE FCONE FCONE FCONE);

      // var U
      maternRaster(
          Axmin, Axres, AxN,
          Aymax, Ayres, AyN,
          varGrid,
          &range[Dparam], &shape[Dparam], &variance[Dparam],
          &anisoRatio[Dparam], &anisoAngleRadians[Dparam],
          &oneI);

      // want varU - covDataGrid %*% t(covDataGrid)
      // Linv %*% covUV is small by big
      // covDataGrid is big by small
      /*
       * B = A = covDataGrid, opB = T, C=varGrid, beta=1 alpha=-1
       * transa = 'N' transB= 'T'
       * M=Ngrid N=Ngrid K=Ny
       * alpha=-1.0
       * A=covDataGrid LDA = Ngrid
       * B=covDataGrid LDB = Ngrid
       * beta = 1.0
       * C = varGrid LDC = Ngrid
       *
       * DGEMM  performs one of the matrix-matrix operations

    C := alpha*op( A )*op( B ) + beta*C,

       */
      F77_NAME(dgemm)(
          "N", "T",
          &Ncell, &Ncell, Ny,
          &minusOneD,
          REAL(covDataGrid), &Ncell,
          REAL(covDataGrid), &Ncell,
          &oneD,
          varGrid, &Ncell
          FCONE FCONE);



      // cholesky
      F77_CALL(dpotrf)("L",
          &Ncell, varGrid,
          &Ncell, &D FCONE);
      if(D != 0){
          error("LAPACK dpotrf failed in maternRasterConditional (conditional covariance), info=%d", D);
      }


      // multiply, want L %*% Z
      //B := alpha*op( A )*B,   or   B := alpha*B*op( A )

      F77_NAME(dtrmm)(
          "R", "L", "N", "N",
          &Ncell, Nsim, &oneD,
          varGrid, &Ncell,
          resultHere, &Ncell
          FCONE FCONE FCONE FCONE);


      // conditional mean
      // covUY %*% varY^(-1) %*% data
      // first L %*% data
      F77_NAME(dtrmm)(
          "R", "L", "N", "N",
          Ny, Nsim, &oneD,
          REAL(varY), Ny,
          ydataHere, Ny
          FCONE FCONE FCONE FCONE);

      // crossprod and add random bit
        F77_NAME(dgemm)(
            "N", "N",
            &Ncell, Nsim, Ny,
            &oneD,
            REAL(covDataGrid), 
            &Ncell,
            ydataHere, Ny,
            &oneD,
            resultHere, &Ncell
            FCONE FCONE);
  } // param loop

//  free(varY);
 // free(varGrid);
//  free(covDataGrid);
  UNPROTECT(2);
}
