#ifndef ACE_CLOSURE_OPS_DEFINED_H
#define ACE_CLOSURE_OPS_DEFINED_H

#include <vector>
#include <Eigen/Dense>
#include <iostream>

namespace ACE{
class Parameters;

class Closure_Ops{
public:
  bool use;
  bool use_env_ops;
  std::vector<Eigen::VectorXcd> ops;
  std::vector<std::pair<int, double> > env_ops_nr;
  double H_scale;
  

  void check_ML(int ML)const;
  
  inline operator bool() const{
    return use;
  }

  void set_ops_from_matrices(const std::vector<Eigen::MatrixXcd> &mats);

  void print_info(std::ostream &ofs=std::cout)const;

  void setup(Parameters &param);
  
  void set_trivial();
  
  inline Closure_Ops(Parameters &param){
    setup(param);
  }
  Closure_Ops();
};

}//namespace
#endif
