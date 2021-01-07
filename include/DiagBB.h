#ifndef DIAG_BB_DEFINED_H
#define DIAG_BB_DEFINED_H

#include "Function.h"
#include "Constants.h"
#include "Coupling_Groups.h"

class DiagBB_Config{
public:
  RealFunctionPtr J;
  double temperature;
  double omega_cutoff;
  size_t Ndiscr;
  bool noSubPS;
  std::vector<double> couplings;
  Coupling_Groups groups;

  virtual int get_dim()const{ return couplings.size(); }

  void set_defaults(){
    omega_cutoff=50;
    Ndiscr=1e6;
    if(couplings.size()<2){
      std::cerr<<"DiagBB:set_defaults: couplings.size()<2!"<<std::endl;
      exit(1);
    }
    noSubPS=false;
  }
  DiagBB_Config(const Coupling_Groups & grp, 
                const Eigen::MatrixXcd &couplings_, 
                RealFunctionPtr J_, double temperature_, bool noSubPS_=false)
    : J(J_), temperature(temperature_), groups(grp) {

//    std::cout<<"NGRPS: "<<grp.Ngrps<<std::endl;

    couplings.resize(grp.Ngrps);
    for(int i=0; i<couplings_.rows(); i++){
      couplings[grp[i]]=couplings_(i,i).real();
    }

     set_defaults();
     noSubPS=noSubPS_;
  }
/*
  DiagBB_Config(){
    set_defaults();
  }
*/
};



class DiagBB: public DiagBB_Config{
public:
  
  double get_beta() {
    if(temperature<1e-10)return 0.;
    return 1./(Constants::kB_in_meV_by_K*temperature);
  }


  class K0_integrand_class: public ComplexFunction{
  public:
    RealFunction *J;
    double beta, dt;
    virtual std::complex<double> f(double w)const{
      if(w*w<1e-20)return 0.;
      return J->f(w)/(w*w) *
std::complex<double>(
Constants::coth( 0.5*beta*Constants::hbar_in_meV_ps*w )*(1.-cos(w*dt)) ,
sin(w*dt));
    }
    K0_integrand_class(RealFunction *J_,double beta_, double dt_)
     : J(J_), beta(beta_), dt(dt_){
    }
  };
  class K0_integrand_class_noSubPS: public ComplexFunction{
  public:
    RealFunction *J;
    double beta, dt;
    virtual std::complex<double> f(double w)const{
      if(w*w<1e-20)return 0.;
      return J->f(w)/(w*w) *
std::complex<double>(
Constants::coth( 0.5*beta*Constants::hbar_in_meV_ps*w )*(1.-cos(w*dt)) ,
sin(w*dt)-w*dt );
    }
    K0_integrand_class_noSubPS(RealFunction *J_,double beta_, double dt_)
     : J(J_), beta(beta_), dt(dt_){
    }
  };
  class Kn_integrand_class: public ComplexFunction{
  public:
    RealFunction *J;
    double beta, dt, tau;
    virtual std::complex<double> f(double w)const{
      if(w*w<1e-20)return 0.;
      return 2.*J->f(w)/(w*w)*(1.-cos(w*dt)) *
std::complex<double>(
Constants::coth( 0.5*beta*Constants::hbar_in_meV_ps*w )*cos(w*tau) ,
-sin(w*tau) );

    }
    Kn_integrand_class(RealFunction *J_,double beta_, double dt_, double tau_)
     : J(J_), beta(beta_), dt(dt_), tau(tau_){
    }
  };


  //obtain K_XX for diagonally coupled bath
  std::complex<double> calculate_K(int n, double dt){
    double beta=get_beta();
    if(n==0){
      if(noSubPS){
        K0_integrand_class_noSubPS integrand(J, beta, dt);
        return integrand.integrate(0, omega_cutoff, Ndiscr);
      }else{
        K0_integrand_class integrand(J, beta, dt);
        return integrand.integrate(0, omega_cutoff, Ndiscr);
      }
    }else{
      Kn_integrand_class integrand(J, beta, dt, n*dt);
      return integrand.integrate(0, omega_cutoff, Ndiscr);
    }
  }

  Eigen::MatrixXcd calculate_expS(int n, double dt){

    size_t dim=get_dim();
    size_t NL=dim*dim;
    Eigen::MatrixXcd expS=Eigen::MatrixXcd::Zero(NL, NL);

    std::complex<double> K=calculate_K(n,dt);
    for(int nu_k=0; nu_k<dim; nu_k++){
      for(int mu_k=0; mu_k<dim; mu_k++){
        for(int nu_l=0; nu_l<dim; nu_l++){
          for(int mu_l=0; mu_l<dim; mu_l++){
            std::complex<double> S=
              (-couplings[nu_l]*couplings[nu_k]
               +couplings[nu_l]*couplings[mu_k])*K
             +(-couplings[mu_k]*couplings[mu_l]
               +couplings[nu_k]*couplings[mu_l])*std::conj(K);

            expS(nu_k*dim+mu_k,nu_l*dim+mu_l)=exp(S);
          }
        }
      }
    }
    if(n==0){
      for(int i=0; i<NL; i++){
        for(int j=0; j<NL; j++){  
          if(i!=j)expS(i,j)=0.;
        }
      }
    }
    return expS;
  }



  DiagBB(const Coupling_Groups & groups, 
         const Eigen::MatrixXcd &couplings, 
         RealFunctionPtr J_, double temperature_, bool noSubPS_=false) 
    : DiagBB_Config(groups, couplings, J_, temperature_, noSubPS_){
  }
  DiagBB(const Eigen::MatrixXcd &couplings, 
         RealFunctionPtr J_, double temperature_, bool noSubPS_=false) 
    : DiagBB_Config(Coupling_Groups(couplings), couplings, J_, temperature_,
        noSubPS_){
  } 
};

#endif
