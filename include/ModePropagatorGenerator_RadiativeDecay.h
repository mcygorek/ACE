#ifndef MODE_PROPAGATOR_GENERATOR_RADIATIVE_DECAY_DEFINED_H
#define MODE_PROPAGATOR_GENERATOR_RADIATIVE_DECAY_DEFINED_H

#include "ModePropagatorGenerator.h"
#include "Parameters.h"
#include "Operators.h"
#include "Operators_Boson.h"

class ModePropagatorGenerator_RadiativeDecay: public ModePropagatorGenerator{
public:

  std::vector<std::pair<double, double> > E_g;  
  int M;
  Eigen::MatrixXcd sysop;

  double get_E(int k)const{ return E_g[k].first; }
  double get_g(int k)const{ return E_g[k].second; }

  virtual int get_N()const{ return 2; }
  virtual int get_N_modes()const{ return E_g.size(); }

  static double rate_from_g(double g, int N, double Ediff){
    return 2.*M_PI*Constants::hbar_in_meV_ps*g*g*(N-1)/Ediff;
  }
  static double g_from_rate(double rate, int N, double Ediff){
    return sqrt(rate/(2.*M_PI*Constants::hbar_in_meV_ps*(N-1)/Ediff));
  }

  static bool compare_abs(const std::pair<double,double> &p1,
                          const std::pair<double,double> &p2){
    return fabs(p1.first)>fabs(p2.first);
  } 

  void setup(int Nmod, int Mdim, double Emin, double Emax, double g_){

    M=Mdim;
    if(Nmod<0){ E_g.clear(); return;}
    if(Nmod==1){
      E_g.resize(1); 
      E_g[0].first=(Emax+Emin)/2.;
      E_g[0].second=g_;
      return;
    }
    
    E_g.resize(Nmod);
    double dE=(Emax-Emin)/(Nmod-1);
    for(int i=0; i<Nmod; i++){
      E_g[i].first=Emin+i*dE;
      E_g[i].second=g_;
    }
#ifdef IF_LEADS_NOSORT 
    sort(E_g.begin(), E_g.end(), compare_abs);
#endif
    
    Operators2x2 op;
    sysop=op.ketbra(0,1);
  }
  void setup_Lorentzian(int Nmod, int Mdim, double Emin, double Emax, double g, double gamma){
    setup(Nmod, Mdim, Emin, Emax, g);
    for(size_t i=0; i<E_g.size(); i++){
      double w=E_g[i].first/Constants::hbar_in_meV_ps;
      E_g[i].second=g/M_PI*gamma/(w*w+gamma*gamma);
    }
  }


  virtual std::vector<Eigen::MatrixXcd> get_env_ops() const{
    std::vector<Eigen::MatrixXcd> mats;
    mats.push_back(Operators_Boson::n(M));
//    mats.push_back(Operators_Boson::n(M)*Operators_Boson::n(M));
    return mats;
  }
  virtual void setup(Parameters &param){
    int N_modes=param.get_as_size_t("RadiativeDecay_N_modes", 0);
    int M=param.get_as_size_t("RadiativeDecay_M", 2);

    double E_min=Constants::hbar_in_meV_ps*param.get_as_double("RadiativeDecay_omega_min", 0.);
           E_min=param.get_as_double("RadiativeDecay_E_min", E_min);

    double E_max=Constants::hbar_in_meV_ps*param.get_as_double("RadiativeDecay_omega_max", 0.);
           E_max=param.get_as_double("RadiativeDecay_E_max", E_max);

    double g=param.get_as_double("RadiativeDecay_g", 0.1);
    if(param.is_specified("RadiativeDecay_rate")){
      if(param.is_specified("RadiativeDecay_g")){
        std::cerr<<"Please do not specify both 'RadiativeDecay_g' and 'RadiativeDecay_rate'!"<<std::endl;
        exit(1);
      }
      if(E_max-E_min<1e-8){
        std::cerr<<"Please only specify 'RadiativeDecay_rate' if also 'RadiativeDecay_E_max' > 'RadiativeDecay_E_min'!"<<std::endl;
        exit(1);
      }
      g=g_from_rate(param.get_as_double("RadiativeDecay_rate"), N_modes, E_max-E_min);
    }

    setup(N_modes, M, E_min, E_max, g);

    if(param.is_specified("RadiativeDecay_SysOp")){
      sysop=param.get_as_operator("RadiativeDecay_SysOp");
    }

    for(size_t i=0; i<E_g.size(); i++){
      std::cout<<E_g[i].first<<" "<<E_g[i].second<<std::endl;
    }
  }

  virtual ModePropagatorPtr getModePropagator(int k)const{
    if(k<0||k>=get_N_modes()){
      std::cerr<<"ModePropagatorGenerator_RadiativeDecay: k<0||k>=get_N_modes()!"<<std::endl; 
      exit(1);
    }

    int sysdim=sysop.rows();
    Eigen::MatrixXcd HB_diag=OuterProduct(
         Eigen::MatrixXcd::Identity(sysdim, sysdim), Operators_Boson::n(M));

    Eigen::MatrixXcd HB_base=
         OuterProduct(sysop, Operators_Boson::adagger(M)) +
         OuterProduct(sysop.adjoint(), Operators_Boson::a(M));

    Eigen::MatrixXcd HB = get_E(k)*HB_diag
           +Constants::hbar_in_meV_ps*get_g(k)*HB_base;

    return ModePropagatorPtr(new ModePropagator(sysdim,get_bath_init(k),HB,get_env_ops()));
  }
  virtual Eigen::MatrixXcd get_bath_init(int k)const{
    return Operators_Boson::vacuum(M);
  }

  ModePropagatorGenerator_RadiativeDecay(int Nmod, int M_, 
                                double Emin, double Emax, double g_){
    setup(Nmod, M_, Emin, Emax, g_);
  }
  ModePropagatorGenerator_RadiativeDecay(Parameters &param){
    setup(param);
  }
};

#endif
