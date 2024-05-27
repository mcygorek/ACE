#ifndef MODE_PROPAGATOR_GENERATOR_QDPHONON_DEFINED_H
#define MODE_PROPAGATOR_GENERATOR_QDPHONON_DEFINED_H

#include "ModePropagatorGenerator.hpp"
#include "Operators_Boson.hpp"
#include "SpectralDensity.hpp"
#include "Equilibrium.hpp"

namespace ACE{

class ModePropagatorGenerator_QDPhonon: public ModePropagatorGenerator{
public:

  int M_max;
  double dw;
  double temperature;
  RealFunctionPtr J;
  bool subtract_polaron_shift;

  virtual std::string name()const{return std::string("QDPhonon");}

  double get_E(int k)const{ return hbar_in_meV_ps*(get_N_modes()-k)*dw; }
  double get_g(int k)const{ 
    double w=get_E(k)/hbar_in_meV_ps;
    return sqrt( J->f(w) * dw );
  }

  virtual int get_N()const{ return 2; }
//  virtual int get_N_modes()const{ return N_modes; }


  virtual void setup(Parameters &param){
    setup_default(param);

    temperature=param.get_as_double("temperature", 4.);
    temperature=param.get_as_double("QDPhonon_temperature", temperature);
   

    double E_max=hbar_in_meV_ps*param.get_as_double("QDPhonon_omega_max", 10.);
           E_max=param.get_as_double("QDPhonon_E_max", E_max);
 
     dw=E_max/hbar_in_meV_ps/get_N_modes();

    M_max=param.get_as_size_t("QDPhonon_M_max", 4);

    J=new SpectralDensity_QD();

    subtract_polaron_shift=param.get_as_bool("QDPhonon_subtract_polaron_shift",
                           param.get_as_bool("subtract_polaron_shift", true));
    
    if(param.is_specified(add_name("print_J"))){
      J->print(param.get_as_string(add_name("print_J")),0., 10, 1000);
    }
    if(param.is_specified(add_name("print_E_g"))){
      std::ofstream ofs(param.get_as_string(add_name("print_E_g")).c_str());
      for(int k=0; k<get_N_modes(); k++){
        ofs<<get_E(k)<<" "<<get_g(k)<<std::endl;
      }
    }
  }

  virtual ModePropagatorPtr getModePropagator(int k)const{
    if(k<0||k>=get_N_modes()){
      std::cerr<<"ModePropagatorGenerator_Leads: k<0||k>=get_N_modes()!"<<std::endl; 
      exit(1);
    }
    const double & hbar=hbar_in_meV_ps;

    Operators op(2);
    Eigen::MatrixXcd HB_diag=otimes(op.id(), Operators_Boson::n(M_max));
    Eigen::MatrixXcd HB_base=otimes(op.ketbra(1,1), 
          Operators_Boson::adagger(M_max)+Operators_Boson::a(M_max));
    Eigen::MatrixXcd HB_shift=otimes(op.ketbra(1,1), 
          Eigen::MatrixXcd::Identity(M_max, M_max));

 
    double E=get_E(k);
    double g=get_g(k);

    Eigen::MatrixXcd HB = E*HB_diag + hbar*g*HB_base;

    if(subtract_polaron_shift && fabs(E)>1e-12){
      HB+=(hbar*g*hbar*g/E)*HB_shift;
    }

    return ModePropagatorPtr(new ModePropagator(2,get_bath_init(k),HB));
  }
  virtual Eigen::MatrixXcd get_bath_init(int k)const{
//    double x=get_E(k)/(Constants::kB_in_meV_by_K*temperature);
    return Boson_Equilibrium(M_max,get_E(k),temperature);
  }

  ModePropagatorGenerator_QDPhonon(){
  }
  ModePropagatorGenerator_QDPhonon(Parameters &param){
    setup(param);
  }
};



}//namespace
#endif
