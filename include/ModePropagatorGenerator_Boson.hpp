#ifndef MODE_PROPAGATOR_GENERATOR_BOSON_DEFINED_H
#define MODE_PROPAGATOR_GENERATOR_BOSON_DEFINED_H

#include <vector>
#include <Eigen/Dense>
#include "ModePropagatorGenerator.hpp"
#include "Function.hpp"

namespace ACE{

class ModePropagatorGenerator_Boson: public ModePropagatorGenerator{
public:

  MPG_Discretization_E_g E_g;
  int M;
  int reduced, reduced_fb;
//  int M_base;
  Eigen::MatrixXcd sysop;
  Parameters gparam;
  bool use_anharmonic; double anharmonic_chi;
  double thermalize;

  bool use_initial_thermal; double E_shift_init, temperature;
  bool use_initial_coherent; std::complex<double> initial_coherent;
  bool use_polaron_shift;

  bool use_pos_phase; double pos_phase_offset, pos_phase_pos, pos_phase_c;

  bool use_env_filter;
  RealFunctionPtr env_filter;

  inline virtual std::string name()const{return std::string("Boson");}

  inline double get_E(int k)const{ return E_g.get_E(k); }
  inline double get_g(int k)const{ return E_g.get_g(k); }
  inline virtual void zero_pad(int N_new){ 
        E_g.zero_pad(N_new); skip_list.resize(N_new, false);}

  inline virtual int get_N()const{ 
    if(sysop.rows()<=2) return 2; 
    else return sysop.rows();
  }
  virtual double k_label(int k)const;

  //print initial boson number per mode
  void print_initial_n(const std::string &fname);
  
  double env_ops_filter(int k)const;
  
  virtual std::vector<Eigen::MatrixXcd> get_env_ops(int k)const;
  
  Eigen::MatrixXcd get_HE_diag(int k)const;
  Eigen::MatrixXcd get_HE(int k)const;
  virtual void setup(Parameters &param);

  /* multiply position-dependent phase factor exp(i k*r) to coupling
     Note: here we need to distinguish between phyiscal wave vectors k 
     and mode indices, which we now call "j". Actually, the indices correspond
     to an energy discretization (1d):
     hbar*omega= offset +  get_E 
     k=omega/c
 
  */
  Eigen::MatrixXcd get_pos_phase_op(int j)const;

  virtual ModePropagatorPtr getModePropagator(int k)const;

  virtual Eigen::MatrixXcd get_bath_init(int k)const;

  inline ModePropagatorGenerator_Boson(Parameters &param){
    setup(param);
  }
};

}//namespace
#endif
