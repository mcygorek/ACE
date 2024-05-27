#ifndef INFLUENCE_FUNCTIONAL_REPEAT_DEFINED_H
#define INFLUENCE_FUNCTIONAL_REPEAT_DEFINED_H

#include <vector>
#include "IF_OD_Abstract.hpp"
#include "MPS.hpp"

namespace ACE{

class InfluenceFunctional_Repeat: public IF_OD_Abstract, public MPS{
public:

  //closures:
  std::vector<Eigen::VectorXcd> c;
  //Environment operators:
  std::vector<std::vector<Eigen::VectorXcd> > env_ops;
  
  void print_debug(std::ostream &ofs=std::cout)const;
  
  //Implementation of IF_OD_Abstract
  virtual const MPS_Matrix & get_a(int n)const;
  
  virtual const Eigen::VectorXcd & get_c(int n)const;
  
  virtual const std::vector<Eigen::VectorXcd> &get_env_ops(int n)const;

  virtual void check_within_limits(int n)const;

  //Genuine IF_OD functions:
  void calculate_closures();

  void read_binary(const std::string &fname);
  
  void write_binary(const std::string &fname)const;

  InfluenceFunctional_Repeat(int N=2);
  
  InfluenceFunctional_Repeat(const std::string &fname);

  inline virtual ~InfluenceFunctional_Repeat(){}
};

}//namespace
#endif
