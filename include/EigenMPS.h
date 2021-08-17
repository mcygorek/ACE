#ifndef ACE_EIGEN_MPS_DEFINED_H
#define ACE_EIGEN_MPS_DEFINED_H

#include "Abstract_EigenMPS.h"


class EigenMPS: public Abstract_EigenMPS{
public: 
  std::vector<EigenMPS_Matrix> a;

  virtual int size()const{ return a.size(); }
  virtual const EigenMPS_Matrix & operator[](int n) const{
    if(n<0||n>=size()){
      std::cerr<<"EigenMPS: accessing element with n<0||n>=size(): "<<n<<"/"<<size()<<"!"<<std::endl;
      exit(1);
    }
    return a[n];
  }

  EigenMPS(int nmax=0) : a(nmax){
  }
  EigenMPS(int nmax, int omax) : a(nmax, EigenMPS_Matrix(omax,1,1)){
    for(int n=0; n<nmax; n++){
      a[n].fill(1.);
    } 
  }
  virtual ~EigenMPS(){}
};

#endif
