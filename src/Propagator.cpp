#include "Propagator.hpp"

namespace ACE{

  Eigen::MatrixXcd Propagator::propagate_noupdate(const Eigen::MatrixXcd &rho){
    int dim=rho.rows();
    Eigen::MatrixXcd tmp=Eigen::MatrixXcd::Zero(dim,dim);
    for(int i=0; i<dim; i++){
      for(int j=0; j<dim; j++){
        for(int k=0; k<dim; k++){
          for(int l=0; l<dim; l++){
            tmp(i,j)+=M(i*dim+j, k*dim+l)*rho(k,l);
          }
        }
      }
    }
    return tmp;
  }

  Eigen::MatrixXcd Propagator::propagate(const Eigen::MatrixXcd &rho, double t, double dt){
    update(t, dt);
    return propagate_noupdate(rho);
  }

  Eigen::MatrixXcd Propagator::get_Htot(double t)const{
    return Eigen::MatrixXcd::Zero(get_dim(), get_dim());
  }
}//namespace
