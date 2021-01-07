#ifndef MODE_PROPAGATOR_GENERATOR_SUPERRADIANCE_DEFINED_H
#define MODE_PROPAGATOR_GENERATOR_SUPERRADIANCE_DEFINED_H

#include "ModePropagatorGenerator.h"
#include "Operators_Boson.h"
#include "SpectralDensity.h"

class ModePropagatorGenerator_Superradiance: public ModePropagatorGenerator{
public:

  std::complex<double> phasefac1, phasefac2;
  std::vector<std::pair<double, double> > E_g;
  int M_max;

  double get_E(int k)const{ return E_g[k].first; }
  double get_g(int k)const{ return E_g[k].second; }

  virtual int get_N()const{ return 4; }
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

  virtual std::vector<Eigen::MatrixXcd> get_env_ops() const{
    std::vector<Eigen::MatrixXcd> mats(1,Operators_Boson::n(M_max));
    return mats;
  }
  void setup(int Nmod, int Mmax, double Emin, double Emax, double g_){
    
    phasefac1=phasefac2=1.;
    M_max=Mmax;
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
  }
  virtual void setup(Parameters &param){
    M_max=param.get_as_size_t("Superradiance_M_max", 2);
    int N_modes=param.get_as_size_t("Superradiance_N_modes", 0);

    double E_min=Constants::hbar_in_meV_ps*param.get_as_double("Superradiance_omega_min", 0.);
           E_min=param.get_as_double("Superradiance_E_min", E_min);

    double E_max=Constants::hbar_in_meV_ps*param.get_as_double("Superradiance_omega_max", 0.);
           E_max=param.get_as_double("Superradiance_E_max", E_max);

    double g=param.get_as_double("Superradiance_g", 0.1);
    if(param.is_specified("Superradiance_rate")){
      if(param.is_specified("Superradiance_g")){
        std::cerr<<"Please do not specify both 'Superradiance_g' and 'Superradiance_rate'!"<<std::endl;
        exit(1);
      }
      if(E_max-E_min<1e-8){
        std::cerr<<"Please only specify 'Superradiance_rate' if also 'Superradiance_E_max' > 'Superradiance_E_min'!"<<std::endl;
        exit(1);
      }
      g=g_from_rate(param.get_as_double("Superradiance_rate"), N_modes, E_max-E_min);
    }

    setup(N_modes, M_max, E_min, E_max, g);
   
    {
      double phase1=param.get_as_double("Superradiance_phase1",0.);
      phasefac1=exp(std::complex<double>(0, phase1*M_PI));
      double phase2=param.get_as_double("Superradiance_phase2",0.);
      phasefac2=exp(std::complex<double>(0, phase2*M_PI));
    }

    for(size_t i=0; i<E_g.size(); i++){
      std::cout<<E_g[i].first<<" "<<E_g[i].second<<std::endl;
    }
  }

  virtual ModePropagatorPtr getModePropagator(int k)const{
    if(k<0||k>=get_N_modes()){
      std::cerr<<"ModePropagatorGenerator_Superradiance: k<0||k>=get_N_modes()!"<<std::endl; 
      exit(1);
    }
    if(get_N()!=4){
      std::cerr<<"ModePropagatorGenerator_Superradiance: So far only implemented for get_N()==4!"<<std::endl;
      exit(1);
    }

    Eigen::MatrixXcd op_id=Eigen::MatrixXcd::Identity(get_N(), get_N());

    Eigen::MatrixXcd s_plus=
                  phasefac1*OuterProduct(Operators2x2::sigma_plus(), 
                                         Eigen::MatrixXcd::Identity(2,2))
                 +phasefac2*OuterProduct(Eigen::MatrixXcd::Identity(2,2),
                                         Operators2x2::sigma_plus());

    Eigen::MatrixXcd s_minus=
       std::conj(phasefac1)*OuterProduct(Operators2x2::sigma_minus(), 
                                         Eigen::MatrixXcd::Identity(2,2))
      +std::conj(phasefac2)*OuterProduct(Eigen::MatrixXcd::Identity(2,2),
                                         Operators2x2::sigma_minus());



    Eigen::MatrixXcd HB_diag=OuterProduct(op_id, Operators_Boson::n(M_max));
    Eigen::MatrixXcd HB_base=
               OuterProduct(s_minus, Operators_Boson::adagger(M_max)) 
              +OuterProduct(s_plus,  Operators_Boson::a(M_max));

    Eigen::MatrixXcd HB = get_E(k)*HB_diag
           +Constants::hbar_in_meV_ps*get_g(k)*HB_base;

    return ModePropagatorPtr(new ModePropagator(get_N(),get_bath_init(k),HB,get_env_ops()));

  }
  virtual Eigen::MatrixXcd get_bath_init(int k)const{
//    double x=get_E(k)/(Constants::kB_in_meV_by_K*temperature);
//    return Operators_Boson::equilibrium(M_max,x);
    return Operators_Boson::vacuum(M_max);
  }

  ModePropagatorGenerator_Superradiance(){
  }
  ModePropagatorGenerator_Superradiance(Parameters &param){
    setup(param);
  }
};




#endif
