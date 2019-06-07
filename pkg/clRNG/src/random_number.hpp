#include <clRNG/mrg31k3p.h>
#include "clRNG.hpp"


// clRNG -> Matrix
void convertclRngMat(clrngMrg31k3pStream* streams, Rcpp::IntegerMatrix result) {
  
  int numWorkItems = result.nrow();
  int Ditem,Delement,Dcis,Dg;
  for(Ditem =0;Ditem < numWorkItems;Ditem++){
    for(Delement=0;Delement < 3;Delement++){
      
      Dcis=0;
      Dg=0;
      result(Ditem, Dcis*6 + Dg*3 + Delement) = streams[Ditem].current.g1[Delement];//0,0; 0,1
      Dg=1;
      result(Ditem,Dcis*6 + Dg*3 + Delement) = streams[Ditem].current.g2[Delement];//0,3; 0,4
      
      Dcis=1;
      Dg=0;
      result(Ditem,Dcis*6 + Dg*3 + Delement) = streams[Ditem].initial.g1[Delement];//0,6; 0,7
      Dg=1;
      result(Ditem,Dcis*6 + Dg*3 + Delement) = streams[Ditem].initial.g2[Delement];//0,9; 0,10
      
      Dcis=2;
      Dg=0;
      result(Ditem,Dcis*6 + Dg*3 + Delement) = streams[Ditem].substream.g1[Delement];//0,12
      Dg=1;
      result(Ditem,Dcis*6 + Dg*3 + Delement) = streams[Ditem].substream.g2[Delement];//0,15
    }
  }
  
}

//matrix ->clRNG streams
void convertMatclRng(Rcpp::IntegerMatrix Sin, clrngMrg31k3pStream* streams){
  
  int Ditem,Delement,Dcis,Dg;
  int numWorkItems = Sin.nrow();
  
  for(Ditem =0;Ditem < numWorkItems;Ditem++){
    for(Delement=0;Delement < 3;Delement++){
      
      Dcis=0;
      Dg=0;
      streams[Ditem].current.g1[Delement] = Sin(Ditem,Dcis*6 + Dg*3 + Delement);
      Dg=1;
      streams[Ditem].current.g2[Delement] = Sin(Ditem,Dcis*6 + Dg*3 + Delement);
      
      Dcis=1;
      Dg=0;
      streams[Ditem].initial.g1[Delement] = Sin(Ditem,Dcis*6 + Dg*3 + Delement);
      Dg=1;
      streams[Ditem].initial.g2[Delement]=Sin(Ditem,Dcis*6 + Dg*3 + Delement);
      
      Dcis=2;
      Dg=0;
      streams[Ditem].substream.g1[Delement]=Sin(Ditem,Dcis*6 + Dg*3 + Delement);
      Dg=1;
      streams[Ditem].substream.g2[Delement] = Sin(Ditem,Dcis*6 + Dg*3 + Delement);
    }
  }
  
  
}