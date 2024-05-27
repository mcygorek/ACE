#ifndef MODE_PROPAGATOR_GENERATOR_BEC_DEFINED_H
#define MODE_PROPAGATOR_GENERATOR_BEC_DEFINED_H

#include "ModePropagatorGenerator.hpp"
#include <vector>
#include <Eigen/Dense>

namespace ACE{

class ModePropagatorGenerator_BEC: public ModePropagatorGenerator{
public:
  MPG_Discretization_E_g E_g;
  int M; //mode dimiension
  int N; //system dimension
  double temperature;

  inline virtual std::string name()const{return std::string("BEC");}
  inline double get_E(int k)const{ return E_g.get_E(k); }
  inline double get_g(int k)const{ return E_g.get_g(k); }
  inline virtual int get_N()const{ return N; }

  virtual void setup(Parameters &param);

  virtual ModePropagatorPtr getModePropagator(int k)const;
  
  virtual Eigen::MatrixXcd get_bath_init(int k)const;

  inline ModePropagatorGenerator_BEC(Parameters &param){
    setup(param);
  }
};

}//namespace
#endif
