#ifndef EQUILIBRIUM_DEFINED_H
#define EQUILIBRIUM_DEFINED_H

#include <Eigen/Core>
#include <iostream>
#include "Constants.h"


Eigen::MatrixXcd Equilibrium(const Eigen::VectorXcd &E, double T){ 
  int dim=E.size();
  if(T<1e-12){  //zero temperature:
    //determine degeneracy:
    int n=1;
    for(int i=1; i<dim; i++){
      if(abs(E(i)-E(0))<1e-12)n++;
    }
    Eigen::MatrixXcd mat=Eigen::MatrixXcd::Zero(dim,dim);
    for(int i=0; i<n; i++){
      mat(i,i)=1./( (double)n);
    }
    return mat;
  }else if(T>=1e6){ // infinite temperature
    Eigen::MatrixXcd mat=1./((double)dim)*Eigen::MatrixXcd::Identity(dim,dim);
    return mat;
  }else{ // finite temperature
    double beta=1./(Constants::kB_in_meV_by_K*T);
    Eigen::MatrixXcd mat=Eigen::MatrixXcd::Zero(dim,dim);
    double norm=0.;
    for(int i=0; i<dim; i++){
      double contrib=exp(-beta*(E(i).real()-E(0).real()));
      mat(i,i)=contrib;
      norm+=contrib;
    }
    for(int i=0; i<dim; i++){
      mat(i,i)/=norm;
    }
    return mat;
  }
}

/*
  Eigen::MatrixXcd Equilibrium(const Eigen::MatrixXcd &H, double T){ 
    int n_max=H.rows();
    if(H.rows()!=H.cols()){
      std::cerr<<"Operators_Boson::equilibrium: H.rows()!=H.cols()!"<<std::endl;
      exit(1);
    }
    if(H.rows()<1){
      std::cerr<<"Operators_Boson::equilibrium: H.rows()<1!"<<std::endl;
      exit(1);
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
*/

#endif
