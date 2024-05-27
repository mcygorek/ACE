#pragma once
#ifndef ACE_OPERATORS_DEFINED_H_
#define ACE_OPERATORS_DEFINED_H_
#include <Eigen/Core>

namespace ACE{

class Operators{
public:
  int N;

  inline Eigen::MatrixXcd ketbra(int i, int j)const{
    Eigen::MatrixXcd mat=Eigen::MatrixXcd::Zero(N,N);
    mat(i,j)=1;
    return mat;
  }

  inline Eigen::MatrixXcd zero()const{ return Eigen::MatrixXcd::Zero(N,N); }

  inline Eigen::MatrixXcd id()const{ return Eigen::MatrixXcd::Identity(N,N); }

  Operators(int N_=0) : N(N_){}
};   


  extern Eigen::MatrixXcd sigma_x();
  extern Eigen::MatrixXcd sigma_y();
  extern Eigen::MatrixXcd sigma_z();
  
  extern Eigen::MatrixXcd sigma_plus();
  extern Eigen::MatrixXcd sigma_minus();
};


#endif
