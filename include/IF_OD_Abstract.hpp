#ifndef IF_OD_ABSTRACT_DEFINED_H
#define IF_OD_ABSTRACT_DEFINED_H

#include "otimes.hpp"
#include "TimeGrid.hpp"
#include "IF_OD_Dictionary.hpp"

/** Used as a wrapper to map more general objects of IF_ODs for 
    the purpose of propagating a density matrix */

namespace ACE{

class IF_OD_Abstract{
public:
  TimeGrid tgrid;
  IF_OD_Dictionary dict;

  virtual const MPS_Matrix & get_a(int n)const=0;
  virtual const Eigen::VectorXcd & get_c(int n)const=0;
  virtual const std::vector<Eigen::VectorXcd> & get_env_ops(int n)const=0;

  //Check if MPO chain is long enough
  virtual void check_within_limits(int n)const=0;

  inline virtual int get_sys_dim()const{
    return dict.get_N();
  }

  virtual ~IF_OD_Abstract(){}
};

}//namespace
#endif
