#ifndef MODE_PROPAGATOR_GENERATOR_SYSTEM_CAVITY_QDPHONON_DEFINED_H
#define MODE_PROPAGATOR_GENERATOR_SYSTEM_CAVITY_QDPHONON_DEFINED_H

#include "ModePropagatorGenerator_QDPhonon.h"
#include "Operators_Boson.h"
#include "SpectralDensity.h"

class ModePropagatorGenerator_system_cavity_QDPhonon: public ModePropagatorGenerator_QDPhonon{
public:

  int n_dim_cavity;

  virtual int get_N()const{ return 2*n_dim_cavity; }


  virtual void setup(Parameters &param){
    n_dim_cavity=param.get_as_size_t("add_system_cavity",0)+1;

    int N_modes=param.get_as_size_t("system_cavity_QDPhonon_N_modes", 0);

    temperature=param.get_as_double("temperature", 4.);
    temperature=param.get_as_double("system_cavity_QDPhonon_temperature", temperature);
   

    double E_max=Constants::hbar_in_meV_ps*param.get_as_double("system_cavity_QDPhonon_omega_max", 10.);
           E_max=param.get_as_double("system_cavity_QDPhonon_E_max", E_max);

    int M_max=param.get_as_size_t("system_cavity_QDPhonon_M_max", 4);

    RealFunctionPtr J_=new SpectralDensity_QD();
    ModePropagatorGenerator_QDPhonon::setup(J_, temperature, N_modes, E_max, M_max);

    subtract_polaron_shift=param.get_as_bool("system_cavity_QDPhonon_subtract_polaron_shift", false);
    if(param.is_specified("system_cavity_QDPhonon_printSD")){
      std::ofstream ofs(param.get_as_string("system_cavity_QDPhonon_printSD").c_str());
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
    const double & hbar=Constants::hbar_in_meV_ps;

    Operators2x2 op;
    Eigen::MatrixXcd cav_id=Eigen::MatrixXcd::Identity(n_dim_cavity,n_dim_cavity);

    Eigen::MatrixXcd HB_diag=OuterProduct(OuterProduct(op.id(),cav_id), Operators_Boson::n(M_max));
    Eigen::MatrixXcd HB_base=OuterProduct(OuterProduct(op.ketbra(1,1),cav_id), 
          Operators_Boson::adagger(M_max)+Operators_Boson::a(M_max));
    Eigen::MatrixXcd HB_shift=OuterProduct(OuterProduct(op.ketbra(1,1),cav_id), 
          Eigen::MatrixXcd::Identity(M_max, M_max));

 
    double E=get_E(k);
    double g=get_g(k);

    Eigen::MatrixXcd HB = E*HB_diag + hbar*g*HB_base;

    if(subtract_polaron_shift && fabs(E)>1e-12){
      HB+=(hbar*g*hbar*g/E)*HB_shift;
    }

    return ModePropagatorPtr(new ModePropagator(get_N(),get_bath_init(k),HB));
  }
  virtual Eigen::MatrixXcd get_bath_init(int k)const{
    double x=get_E(k)/(Constants::kB_in_meV_by_K*temperature);
    return Operators_Boson::equilibrium(M_max,x);
  }

  ModePropagatorGenerator_system_cavity_QDPhonon(){
  }
  ModePropagatorGenerator_system_cavity_QDPhonon(Parameters &param){
    setup(param);
  }
};




#endif
