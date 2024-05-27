#ifndef REP_GIF_DEFINED_H
#define REP_GIF_DEFINED_H

#include <vector>
#include <Eigen/Dense>
#include "MPS_Matrix.hpp"
#include "IF_OD_Abstract.hpp"
#include "otimes.hpp"
#include "ModePropagator.hpp"
#include "Compress_Trafo_At.hpp"

namespace ACE{

class InfluenceFunctional_OD;

class Rep_GIF: public IF_OD_Abstract{
public:

  MPS_Matrix M;
  Eigen::VectorXcd init;
  std::vector<Eigen::VectorXcd> env_ops;
 
  virtual const MPS_Matrix &get_a(int n)const{
    return M;
  }
  virtual const Eigen::VectorXcd &get_c(int n)const{
    return env_ops[0];
  }
  virtual const std::vector<Eigen::VectorXcd> &get_env_ops(int n)const{
    return env_ops;
  }

  virtual void check_within_limits(int n)const{ 
    return;
  }

  void set_default(int N);
  
  void expand_ops(const ModePropagator &mprop);

  void apply_compress_trafo(Compress_Trafo_At &cta, bool low_high_low);

  void regularize();
 
  inline void read_binary(const std::string &fname){
  }
  inline void write_binary(const std::string &fname){
  }
   
  inline Rep_GIF(int N=2){
    set_default(N);
  }
};

}//namespace
#endif
