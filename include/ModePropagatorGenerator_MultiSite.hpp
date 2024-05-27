#ifndef MODE_PROPAGATOR_GENERATOR_MULTISITE_DEFINED_H
#define MODE_PROPAGATOR_GENERATOR_MULTISITE_DEFINED_H

#include "ModePropagatorGenerator.hpp"

namespace ACE{
class ModePropagatorGenerator_MultiSite: public ModePropagatorGenerator{
public:

  MPG_Discretization_E_g E_g;
  int M;
  int Nintermediate;

  std::vector<double> tpos;

  int reduced, reduced_fb;
  Parameters gparam;

  bool use_initial_thermal; double EFermi, temperature;
  bool use_initial_coherent; std::complex<double> initial_coherent;


  bool use_IP;
  Eigen::MatrixXcd H_IP;

  inline virtual std::string name()const{return std::string("MultiSite");}

  inline int get_N_sites()const{ return tpos.size();}
  inline double get_E(int k)const{ return E_g.get_E(k); }
  inline double get_omega(int k)const{ return E_g.get_E(k)/hbar_in_meV_ps; }
  inline double get_g(int k)const{ return E_g.get_g(k); }

  inline virtual int get_N()const{ 
    if(get_N_sites()<=2) return 2; 
    else return get_N_sites();
  }

  //odd-even 
  inline double phase_sign(int k)const{
    if(k%2==0)return 1.;
    else return -1.;
  }

  //print initial boson number per mode
  void print_initial_n(const std::string &fname);
  
  virtual std::vector<Eigen::MatrixXcd> get_env_ops(int k) const;
  
  virtual void setup(Parameters &param);

  virtual ModePropagatorPtr getModePropagator(int k)const;

  virtual Eigen::MatrixXcd get_bath_init(int k)const;

  inline ModePropagatorGenerator_MultiSite(Parameters &param){
    setup(param);
  }
};
}//namespace
#endif
