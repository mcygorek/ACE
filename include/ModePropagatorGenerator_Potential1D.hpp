#ifndef MODE_PROPAGATOR_GENERATOR_POTENTIAL1D_DEFINED_H
#define MODE_PROPAGATOR_GENERATOR_POTENTIAL1D_DEFINED_H

#include "ModePropagatorGenerator.hpp"
#include "Potential1D.hpp"

namespace ACE{
class Parameters;

class ModePropagatorGenerator_Potential1D: public ModePropagatorGenerator{
public:

  MPG_Discretization_E_g E_g;

  Potential1D p1d;
  int M;
  int reduced, reduced_fb;
//  int M_base;
  Eigen::MatrixXcd sysop;
  Parameters gparam;

  double temperature;
//  bool use_polaron_shift; 
  bool use_renorm_x;


  inline virtual std::string name()const{return std::string("Potential1D");}

  inline double get_E(int k)const{ return E_g.get_E(k); }
  inline double get_g(int k)const{ return E_g.get_g(k); }

  inline virtual int get_N()const{ 
    if(sysop.rows()<=2) return 2; 
    else return sysop.rows();
  }

  //print initial boson number per mode
  void print_initial_n(const std::string &fname)const;
  
  virtual EnvironmentOperators get_env_ops(int k)const;
  
  virtual void setup(Parameters &param);

  virtual ModePropagatorPtr get_ModePropagator(int k)const;

  virtual Eigen::MatrixXcd get_bath_init(int k)const;

  inline ModePropagatorGenerator_Potential1D(Parameters &param){
    setup(param);
  }
};

}//namespace
#endif
