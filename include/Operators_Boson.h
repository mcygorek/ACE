#ifndef OPERATORS_BOSON_DEFINED_H
#define OPERATORS_BOSON_DEFINED_H

#include <Eigen/Core>
#include <iostream>

namespace Operators_Boson{

  Eigen::MatrixXcd a(int n_max){
    Eigen::MatrixXcd result=Eigen::MatrixXcd::Zero(n_max, n_max);
    for(int i=0; i<n_max-1; i++){
      result(i, i+1)=sqrt(i+1);
    }
    return result;
  }

  Eigen::MatrixXcd adagger(int n_max){
    Eigen::MatrixXcd result=Eigen::MatrixXcd::Zero(n_max, n_max);
    for(int i=0; i<n_max-1; i++){
      result(i+1, i)=sqrt(i+1);
    }
    return result;
  }

  Eigen::MatrixXcd n(int n_max){
    Eigen::MatrixXcd result=Eigen::MatrixXcd::Zero(n_max, n_max);
    for(int i=0; i<n_max; i++){
      result(i, i)=i;
    }
    return result;
  } 

  Eigen::MatrixXcd equilibrium(int n_max, double x){ //x=(E-mu)/(kB T)
    if(x<1e-8){
      return 1./(double)n_max*Eigen::MatrixXcd::Identity(n_max, n_max);
    }
    Eigen::VectorXd v(n_max);
    double norm=0;
    for(int i=0; i<n_max; i++){
      double e=exp(-x*i);
      norm+=e;
      v(i)=e;
    }
    for(int i=0; i<n_max; i++){
      v(i)/=norm;
    }
    return v.asDiagonal();
  }
  Eigen::MatrixXcd vacuum(int n_max){
    Eigen::MatrixXcd mat=Eigen::MatrixXcd::Zero(n_max, n_max);
    mat(0,0)=1;
    return mat;
  }
};

#endif
