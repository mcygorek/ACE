#ifndef MODE_PROPAGATOR_GENERATOR_QDPHONON_DEFINED_H
#define MODE_PROPAGATOR_GENERATOR_QDPHONON_DEFINED_H

#include "ModePropagatorGenerator.h"
#include "Operators_Boson.h"
#include "SpectralDensity.h"

class ModePropagatorGenerator_QDPhonon: public ModePropagatorGenerator{
public:

  int M_max;
  int N_modes;
  double dw;
  double temperature;
  RealFunctionPtr J;
  bool subtract_polaron_shift;

  double get_E(int k)const{ return Constants::hbar_in_meV_ps*(N_modes-k)*dw; }
  double get_g(int k)const{ 
    double w=get_E(k)/Constants::hbar_in_meV_ps;
    return sqrt( J->f(w) * dw );
  }

  virtual int get_N()const{ return 2; }
  virtual int get_N_modes()const{ return N_modes; }


  void setup(RealFunctionPtr J_, double T, int Nmod, double Emax, int M_max_){
    J=J_;
    temperature=T;
    N_modes=Nmod;
    dw=Emax/Nmod/Constants::hbar_in_meV_ps;
    M_max=M_max_;
//std::cout<<"E_max: "<<Emax<<" "<<get_N_modes()*dw*Constants::hbar_in_meV_ps<<" "<<get_E(0)<<std::endl;

    if(Nmod<1){ std::cerr<<"QDPhonon: N_modes<1!"<<std::endl; exit(1); }
  }
  virtual void setup(Parameters &param){
    int N_modes=param.get_as_size_t("QDPhonon_N_modes", 0);

    temperature=param.get_as_double("temperature", 4.);
    temperature=param.get_as_double("QDPhonon_temperature", temperature);
   

    double E_max=Constants::hbar_in_meV_ps*param.get_as_double("QDPhonon_omega_max", 10.);
           E_max=param.get_as_double("QDPhonon_E_max", E_max);

    int M_max=param.get_as_size_t("QDPhonon_M_max", 4);

    RealFunctionPtr J_=new SpectralDensity_QD();
    setup(J_, temperature, N_modes, E_max, M_max);

    subtract_polaron_shift=param.get_as_bool("QDPhonon_subtract_polaron_shift",
                           param.get_as_bool("subtract_polaron_shift", true));
    
    if(param.is_specified("QDPhonon_printSD")){
      std::ofstream ofs(param.get_as_string("QDPhonon_printSD").c_str());
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
    Eigen::MatrixXcd HB_diag=OuterProduct(op.id(), Operators_Boson::n(M_max));
    Eigen::MatrixXcd HB_base=OuterProduct(op.ketbra(1,1), 
          Operators_Boson::adagger(M_max)+Operators_Boson::a(M_max));
    Eigen::MatrixXcd HB_shift=OuterProduct(op.ketbra(1,1), 
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
    double x=get_E(k)/(Constants::kB_in_meV_by_K*temperature);
    return Operators_Boson::equilibrium(M_max,x);
  }

  ModePropagatorGenerator_QDPhonon(){
  }
  ModePropagatorGenerator_QDPhonon(Parameters &param){
    setup(param);
  }
};




#endif
