#ifndef MODE_PROPAGATOR_GENERATOR_FERMION_DEFINED_H
#define MODE_PROPAGATOR_GENERATOR_FERMION_DEFINED_H

#include "ModePropagatorGenerator.hpp"

namespace ACE{
class Parameters;

class ModePropagatorGenerator_Fermion: public ModePropagatorGenerator{
public:

  MPG_Discretization_E_g E_g;
  int continuum_subdiv_N;
  ModePropagator::low_pass_struct low_pass;

  double EFermi, temperature;
  double thermalize;
  Eigen::MatrixXcd sysop;
  Parameters gparam;

  bool no_env_ops;

  inline virtual std::string name()const{return std::string("Fermion");}

  inline double get_E(int k)const{ return E_g.get_E(k); }
  inline double get_dE(int k)const{ return E_g.get_dE(k); }
  inline double get_g(int k)const{ return E_g.get_g(k); }
  inline virtual void zero_pad(int N_new){ 
        E_g.zero_pad(N_new); skip_list.resize(N_new, false);}

  inline virtual int get_N()const{ 
    if(sysop.rows()<=2) return 2; 
    else return sysop.rows();
  }
  virtual double k_label(int k)const;

  virtual std::vector<Eigen::MatrixXcd> get_env_ops(int k) const;

  virtual void setup(Parameters &param);

  virtual ModePropagatorPtr getModePropagator(int k)const;
  
  double get_n_eq(int k)const; //particle number in therm. equil. 
  virtual Eigen::MatrixXcd get_bath_init(int k)const;

  inline ModePropagatorGenerator_Fermion(Parameters &param){
    setup(param);
  }
  inline virtual ~ModePropagatorGenerator_Fermion(){}
};

}//namespace
#endif
