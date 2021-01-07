#ifndef IF_OD_ABSTRACT_DEFINED_H
#define IF_OD_ABSTRACT_DEFINED_H

#include "OuterProduct.h"
#include "IF_TimeGrid.h"
#include "IF_OD_Dictionary.h"

/** Used as a wrapper to map more general objects of IF_ODs for 
    the purpose of propagating a density matrix */


class IF_OD_Abstract{
public:
  IF_TimeGrid tgrid;
  IF_OD_Dictionary dict;

  virtual const MPS_Matrix & get_a(int n)const=0;
  virtual const Eigen::VectorXcd & get_c(int n)const=0;
  virtual const std::vector<Eigen::VectorXcd> & get_env_ops(int n)const=0;

  //Check if MPO chain is long enough
  virtual void check_within_limits(int n)const=0;


  virtual ~IF_OD_Abstract(){}
};

#endif
