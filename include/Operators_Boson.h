#ifndef OPERATORS_BOSON_DEFINED_H
#define OPERATORS_BOSON_DEFINED_H

#include <Eigen/Core>
#include <iostream>
#include "Constants.h"
#include "Operators.h"

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
  Eigen::MatrixXcd id(int N){ 
    return Eigen::MatrixXcd::Identity(N,N); 
  }

  Eigen::MatrixXcd vacuum(int n_max){
    Eigen::MatrixXcd mat=Eigen::MatrixXcd::Zero(n_max, n_max);
    mat(0,0)=1;
    return mat;
  }
  Eigen::MatrixXcd equilibrium(int n_max, double dE, double T){ //x=(E-mu)
    if(T<1e-12){  //effectively: step function
      if(dE<0){
        return 1./(double)n_max*Eigen::MatrixXcd::Identity(n_max, n_max);
      }else{
        return vacuum(n_max);
      }
    }   

    double x=dE/(Constants::kB_in_meV_by_K*T);
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

  // Equilibrium w.r.t. to Hamiltonian matrix
  Eigen::MatrixXcd equilibrium(Eigen::MatrixXcd H, double T, double E_shift=0.){ 
    int n_max=H.rows();
    if(H.rows()!=H.cols()){
      std::cerr<<"Operators_Boson::equilibrium: H.rows()!=H.cols()!"<<std::endl;
      exit(1);
    }
    if(H.rows()<1){
      std::cerr<<"Operators_Boson::equilibrium: H.rows()<1!"<<std::endl;
      exit(1);
    }

    if(T<1e-20){
      Eigen::MatrixXcd ret=Eigen::MatrixXcd::Zero(H.rows(),H.rows());
      ret(0,0)=1.;
      return ret;
    }

    for(int i=0; i<H.rows(); i++){
      H(i,i)+=E_shift;
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver(H);
    Eigen::MatrixXcd V=solver.eigenvectors();

    int nrltz=0;
    for(int i=0; i<n_max; i++){
      if(solver.eigenvalues()(i)<0.)nrltz++;
    }
    if(nrltz>0){
      std::cerr<<"Operators_Boson::equilibrium: "<<nrltz<<" of "<<n_max<<" eigenvalues smaller than 0!"<<std::endl;
      std::cerr<<solver.eigenvalues()<<std::endl;
      exit(1);
    }

    if(T<1e-12){  //effectively: step function
      return vacuum(n_max);
    }   

    double beta=1./(Constants::kB_in_meV_by_K*T);

    Eigen::VectorXd v(n_max);
    double norm=0;
    for(int i=0; i<n_max; i++){
      double e=exp(-beta*solver.eigenvalues()(i));
      norm+=e;
      v(i)=e;
    }
    for(int i=0; i<n_max; i++){
      v(i)/=norm;
    }
 
    Eigen::MatrixXcd ret=Eigen::MatrixXcd::Zero(n_max,n_max);
    for(int i=0; i<n_max; i++){
       ret += v(i) * V.col(i)*( V.col(i).adjoint());
    }
    return ret;
  }
};

//
namespace Operators_Boson_Offset{

  Eigen::MatrixXcd a(int base, int M){
    Eigen::MatrixXcd result=Eigen::MatrixXcd::Zero(M, M);
    for(int i=0; i<M-1; i++){
      result(i, i+1)=sqrt(base+i+1);
    }
    return result;
  }

  Eigen::MatrixXcd adagger(int base, int M){
    Eigen::MatrixXcd result=Eigen::MatrixXcd::Zero(M,M);
    for(int i=0; i<M-1; i++){
      result(i+1, i)=sqrt(base+i+1);
    }
    return result;
  }

  Eigen::MatrixXcd n(int base, int M){
    Eigen::MatrixXcd result=Eigen::MatrixXcd::Zero(M, M);
    for(int i=0; i<M; i++){
      result(i, i)=base+i;
    }
    return result;
  } 

  Eigen::MatrixXcd lowest(int M){
    Eigen::MatrixXcd mat=Eigen::MatrixXcd::Zero(M, M);
    mat(0,0)=1;
    return mat;
  }
};

#endif
