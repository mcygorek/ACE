#ifndef MODE_PROPAGATOR_GENERATOR_SINGLEMODE_DEFINED_H
#define MODE_PROPAGATOR_GENERATOR_SINGLEMODE_DEFINED_H

#include "ModePropagatorGenerator.hpp"
#include <vector>
#include <Eigen/Dense>

namespace ACE{

class ModePropagatorGenerator_SingleMode: public ModePropagatorGenerator{
public:

  Eigen::MatrixXcd HE;
  Eigen::MatrixXcd rho_init;

  std::vector<Eigen::MatrixXcd> envops;

  inline virtual std::string name()const{return std::string("SingleMode");}

  inline virtual int get_N()const{ return HE.rows()/rho_init.rows(); }

  void check_validity()const;

  void setup(const Eigen::MatrixXcd & HE_,  const Eigen::MatrixXcd & init_);
  
  void setup(const std::vector<Eigen::MatrixXcd> & ops);
  
  void setup(const std::vector<std::string> & str);

  inline virtual std::vector<Eigen::MatrixXcd> get_env_ops(int k) const{
    return envops;
  }
  inline virtual Eigen::MatrixXcd get_bath_init(int k)const{
    return rho_init;
  }
  virtual ModePropagatorPtr getModePropagator(int k)const;

  inline ModePropagatorGenerator_SingleMode(const Eigen::MatrixXcd & HE_,  const Eigen::MatrixXcd & init_){
    setup(HE_, init_);
  }
  inline ModePropagatorGenerator_SingleMode(const std::vector<Eigen::MatrixXcd> &ops){
    setup(ops);
  }
  inline ModePropagatorGenerator_SingleMode(const std::vector<std::string> &svec){
    setup(svec);
  }
};

}//namespace
#endif
