#pragma once
#ifndef ACE_OTIMES_DEFINED_H
#define ACE_OTIMES_DEFINED_H

#include <unsupported/Eigen/KroneckerProduct>

namespace ACE{

inline Eigen::MatrixXcd otimes(const Eigen::MatrixXcd &o1, 
                              const Eigen::MatrixXcd &o2){

  return Eigen::kroneckerProduct(o1,o2);
}
inline Eigen::MatrixXcd otimes(const Eigen::MatrixXcd &o1, 
                              const Eigen::MatrixXcd &o2,
                              const Eigen::MatrixXcd &o3){

  return otimes(otimes(o1,o2),o3);
}

inline Eigen::MatrixXd otimes_real(const Eigen::MatrixXd &o1, 
                              const Eigen::MatrixXd &o2){

  return Eigen::kroneckerProduct(o1,o2);
}

inline Eigen::MatrixXcd ExpandSingleOp(int dim_before, const Eigen::MatrixXcd &m, int dim_after){
  if(dim_before<=1 && dim_after<=1)return m;
  if(dim_before<=1)return otimes(m, Eigen::MatrixXcd::Identity(dim_after,dim_after));
  if(dim_after<=1)return otimes(Eigen::MatrixXcd::Identity(dim_before,dim_before),m);
  return otimes(Eigen::MatrixXcd::Identity(dim_before,dim_before), m,
                      Eigen::MatrixXcd::Identity(dim_after,dim_after));
}


inline Eigen::VectorXcd Vector_otimes(const Eigen::VectorXcd &o1, 
                               const Eigen::VectorXcd &o2){

  return Eigen::kroneckerProduct(o1,o2);
}

}

#endif
