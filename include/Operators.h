#ifndef OPERATORS2X2_DEFINED_H_
#define OPERATORS2X2_DEFINED_H_
#include <Eigen/Core>


class Operators{
public:
  int N;

  Eigen::MatrixXcd ketbra(int i, int j){
    Eigen::MatrixXcd mat=Eigen::MatrixXcd::Zero(N,N);
    mat(i,j)=1;
    return mat;
  }

  Eigen::MatrixXcd zero(){ return Eigen::MatrixXcd::Zero(N,N); }

  Eigen::MatrixXcd id(){ return Eigen::MatrixXcd::Identity(N,N); }

  Operators(int N_=0) : N(N_){}
};   

class Operators2x2{
public:
  static Eigen::MatrixXcd ketbra(int i, int j){
    Eigen::MatrixXcd mat=Eigen::MatrixXcd::Zero(2,2);
    mat(i,j)=1;
    return mat;
  }
  static Eigen::MatrixXcd zero(){ return Eigen::MatrixXcd::Zero(2,2); }
  static Eigen::MatrixXcd id(){ return Eigen::MatrixXcd::Identity(2,2); }

  static Eigen::MatrixXcd sigma_x(){
    Eigen::MatrixXcd mat=Eigen::MatrixXcd::Zero(2,2);
    mat(0,1)=1.;
    mat(1,0)=1.;
    return mat;
  }
  static Eigen::MatrixXcd sigma_y(){
    Eigen::MatrixXcd mat=Eigen::MatrixXcd::Zero(2,2);
    mat(0,1)=std::complex<double>(0.,1.);
    mat(1,0)=std::complex<double>(0.,-1.);
    return mat;
  }
  static Eigen::MatrixXcd sigma_z(){
    Eigen::MatrixXcd mat=Eigen::MatrixXcd::Zero(2,2);
    mat(0,0)=-1.;
    mat(1,1)=1.;
    return mat;
  }
  static Eigen::MatrixXcd sigma_plus(){
    Eigen::MatrixXcd mat=Eigen::MatrixXcd::Zero(2,2);
    mat(1,0)=1;
    return mat;
  }
  static Eigen::MatrixXcd sigma_minus(){
    Eigen::MatrixXcd mat=Eigen::MatrixXcd::Zero(2,2);
    mat(0,1)=1;
    return mat;
  }
};


Eigen::MatrixXcd otimes(const Eigen::MatrixXcd &A, const Eigen::MatrixXcd &B){
  Eigen::MatrixXcd C=Eigen::MatrixXcd::Zero(A.rows()*B.rows(), A.cols()*B.cols());
  for(int i=0; i<A.rows(); i++){
    for(int j=0; j<A.cols(); j++){
      for(int k=0; k<B.rows(); k++){
        for(int l=0; l<B.cols(); l++){
          C(i*B.rows()+k, j*B.cols()+l)=A(i,j)*B(k,l);
        }
      }
    }
  }
  return C;
}

#include "Operators_DotCavity.h"
#include "Operators_Multi2lvl.h"
#include "OuterProduct.h"

#endif
