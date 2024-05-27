#ifndef MODE_PROPAGATOR_GENERATOR_LASING_DEFINED_H
#define MODE_PROPAGATOR_GENERATOR_LASING_DEFINED_H

#include "ModePropagatorGenerator.hpp"
#include <vector>
#include <Eigen/Dense>

namespace ACE{

class ModePropagatorGenerator_Lasing: public ModePropagatorGenerator{
public:

  MPG_Discretization_E_g E_g;
  int N_system;
  double g, g_prime;
  double Gamma_up, Gamma_down;

  inline virtual std::string name()const{return std::string("Lasing");}

  inline virtual int get_N()const{ return N_system; }

  virtual std::vector<Eigen::MatrixXcd> get_env_ops(int k)const;

  virtual void setup(Parameters &param);

  virtual ModePropagatorPtr getModePropagator(int k)const;
  
  virtual Eigen::MatrixXcd get_bath_init(int k)const;

  inline ModePropagatorGenerator_Lasing(Parameters &param){
    setup(param);
  }
};

}//namespace
#endif
