#ifndef MODE_PROPAGATOR_GENERATOR_LEADS_DEFINED_H
#define MODE_PROPAGATOR_GENERATOR_LEADS_DEFINED_H

#include "ModePropagatorGenerator.h"
#include "Parameters.h"
#include "Operators.h"

class ModePropagatorGenerator_Leads: public ModePropagatorGenerator{
public:

  std::vector<std::pair<double, double> > E_g;  
  double EFermi, temperature;

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

  void setup(double EF_, double T, int Nmod, double Emin, double Emax, double g_){
    EFermi=EF_;
    temperature=T;

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

  virtual std::vector<Eigen::MatrixXcd> get_env_ops() const{
    Operators2x2 op; 
    std::vector<Eigen::MatrixXcd> mats;
    mats.push_back(op.ketbra(1,1));
    mats.push_back(op.ketbra(0,1));
    mats.push_back(op.ketbra(1,0));
    return mats;
  }

  virtual void setup(Parameters &param){
    int N_modes=param.get_as_size_t("Leads_N_modes", 0);

    EFermi=param.get_as_double("EFermi", -1e6);
    EFermi=param.get_as_double("Leads_EFermi", EFermi);

    temperature=param.get_as_double("temperature", 4.);
    temperature=param.get_as_double("Leads_temperature", temperature);
   
    double Leads_E_min=Constants::hbar_in_meV_ps*param.get_as_double("Leads_omega_min", 0.);
           Leads_E_min=param.get_as_double("Leads_E_min", Leads_E_min);

    double Leads_E_max=Constants::hbar_in_meV_ps*param.get_as_double("Leads_omega_max", 0.);
           Leads_E_max=param.get_as_double("Leads_E_max", Leads_E_max);

    double Leads_g=param.get_as_double("Leads_g", 0.1);
    if(param.is_specified("Leads_rate")){
      if(param.is_specified("Leads_g")){
        std::cerr<<"Please do not specify both 'Leads_g' and 'Leads_rate'!"<<std::endl;
        exit(1);
      }
      if(Leads_E_max-Leads_E_min<1e-8){
        std::cerr<<"Please only specify 'Leads_rate' if also 'Leads_E_max' > 'Leads_E_min'!"<<std::endl;
        exit(1);
      }
      Leads_g=g_from_rate(param.get_as_double("Leads_rate"), N_modes, Leads_E_max-Leads_E_min);
    }


    setup(EFermi, temperature, N_modes, Leads_E_min, Leads_E_max, Leads_g);

    for(size_t i=0; i<E_g.size(); i++){
      std::cout<<E_g[i].first<<" "<<E_g[i].second<<std::endl;
    }
  }

  virtual ModePropagatorPtr getModePropagator(int k)const{
    if(k<0||k>=get_N_modes()){
      std::cerr<<"ModePropagatorGenerator_Leads: k<0||k>=get_N_modes()!"<<std::endl; 
      exit(1);
    }

    Operators2x2 op;
    Eigen::MatrixXcd HB_diag=OuterProduct(op.id(), op.ketbra(1,1));
    Eigen::MatrixXcd HB_base=OuterProduct(op.ketbra(0,1), op.ketbra(1,0)) + \
                             OuterProduct(op.ketbra(1,0), op.ketbra(0,1));
//    Eigen::MatrixXcd HB_base= \
OuterProduct(op.ketbra(0,1), std::complex<double>(0.,-1.)*op.ketbra(1,0)) + \
OuterProduct(op.ketbra(1,0), std::complex<double>(0., 1.)*op.ketbra(0,1));

//std::complex<double> ephi=exp(std::complex<double>(0,M_PI/2));
//    Eigen::MatrixXcd HB_base=  \
OuterProduct(op.ketbra(0,1), std::conj(ephi)*op.ketbra(1,0)) + \
OuterProduct(op.ketbra(1,0),           ephi *op.ketbra(0,1));

    Eigen::MatrixXcd HB = get_E(k)*HB_diag
           +Constants::hbar_in_meV_ps*get_g(k)*HB_base;

    return ModePropagatorPtr(new ModePropagator(2,get_bath_init(k),HB,get_env_ops()));
  }
  virtual Eigen::MatrixXcd get_bath_init(int k)const{
    if(temperature<1e-6){
      if(get_E(k)>EFermi)return Operators2x2::ketbra(0,0);
      else return Operators2x2::ketbra(1,1);
    }
    double x=get_E(k)/(Constants::kB_in_meV_by_K*temperature);
    double f=Constants::fermi(x);
    return (1.-f)*Operators2x2::ketbra(0,0)+f*Operators2x2::ketbra(1,1);
  }

  ModePropagatorGenerator_Leads(double EF_, double T, int Nmod, 
                                double Emin, double Emax, double g_){
    setup(EF_, T, Nmod, Emin, Emax, g_);
  }
  ModePropagatorGenerator_Leads(Parameters &param){
    setup(param);
  }
};

#endif
