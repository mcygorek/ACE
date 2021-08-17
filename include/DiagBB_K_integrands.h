#ifndef ACE_DIAGBB_K_INTEGRANDS
#define ACE_DIAGBB_K_INTEGRANDS

 
// Conventional integrands:

  class K0_integrand_class: public ComplexFunction{
  public:
    RealFunction *J;
    double beta, E_shift, dt;
    virtual std::complex<double> f(double w)const{
      if(w*w<1e-20)return 0.;
      return J->f(w)/(w*w)*
         std::complex<double>(get_coth(beta,E_shift,w)*(1.-cos(w*dt)) , sin(w*dt));
    }
    K0_integrand_class(RealFunction *J_,double beta_, double E_shift_, double dt_)
     : J(J_), beta(beta_), E_shift(E_shift_), dt(dt_){
    }
  };
  class K0_integrand_class_noSubPS: public ComplexFunction{
  public:
    RealFunction *J;
    double beta, E_shift, dt;
    virtual std::complex<double> f(double w)const{
      if(w*w<1e-20)return 0.;
      return J->f(w)/(w*w) *
std::complex<double>(get_coth(beta,E_shift,w)*(1.-cos(w*dt)) , sin(w*dt)-w*dt );
    }
    K0_integrand_class_noSubPS(RealFunction *J_,double beta_, double E_shift_, double dt_)
     : J(J_), beta(beta_), E_shift(E_shift_), dt(dt_){
    }
  };
  class Kn_integrand_class: public ComplexFunction{
  public:
    RealFunction *J;
    double beta, E_shift, dt, tau;
    virtual std::complex<double> f(double w)const{
      if(w*w<1e-20)return 0.;
      return 2.*J->f(w)/(w*w)*(1.-cos(w*dt)) *
std::complex<double>(get_coth(beta,E_shift,w)*cos(w*tau) , -sin(w*tau) );

    }
    Kn_integrand_class(RealFunction *J_,double beta_, double E_shift_, double dt_, double tau_)
     : J(J_), beta(beta_), E_shift(E_shift_), dt(dt_), tau(tau_){
    }
  };


// new ones: integrate cos and sin parts separately:

  class K0_const_integrand_class: public ComplexFunction{
  public:
    RealFunction *J;
    double beta, E_shift;
    virtual std::complex<double> f(double w)const{
      if(w*w<1e-20)return 0.;
      return J->f(w)/(w*w)*get_coth(beta,E_shift,w);
    }
    K0_const_integrand_class(RealFunction *J_,double beta_, double E_shift_, double dt_)
     : J(J_), beta(beta_), E_shift(E_shift_) {
    }
  };
  class K0_const_integrand_class_noSubPS: public ComplexFunction{
  public:
    RealFunction *J;
    double beta, E_shift, dt;
    virtual std::complex<double> f(double w)const{
      if(w*w<1e-20)return 0.;
      return J->f(w)/(w*w)*std::complex<double>(get_coth(beta,E_shift,w), -w*dt );
    }
    K0_const_integrand_class_noSubPS(RealFunction *J_,double beta_, double E_shift_, double dt_)
     : J(J_), beta(beta_), E_shift(E_shift_), dt(dt_){
    }
  };
   class K0_posdt_integrand_class: public ComplexFunction{
  public:
    RealFunction *J;
    double beta, E_shift;
    virtual std::complex<double> f(double w)const{
      if(w*w<1e-20)return 0.;
      return J->f(w)/(w*w)* 1./2.*( -get_coth(beta,E_shift,w) + 1.);
    }
    K0_posdt_integrand_class(RealFunction *J_,double beta_, double E_shift_)
     : J(J_), beta(beta_), E_shift(E_shift_){
    }
  };
  class K0_negdt_integrand_class: public ComplexFunction{
  public:
    RealFunction *J;
    double beta, E_shift;
    virtual std::complex<double> f(double w)const{
      if(w*w<1e-20)return 0.;
      return J->f(w)/(w*w)* 1./2.*( -get_coth(beta,E_shift,w) - 1.);
    }
    K0_negdt_integrand_class(RealFunction *J_,double beta_, double E_shift_)
     : J(J_), beta(beta_), E_shift(E_shift_){
    }
  };


  class Kn_posdt_integrand_class: public ComplexFunction{
  public:
    RealFunction *J;
    double beta, E_shift, dt;
    virtual std::complex<double> f(double w)const{
      if(w*w<1e-20)return 0.;
      return 2.*J->f(w)/(w*w)*(1.-cos(w*dt))*1./2.*(get_coth(beta,E_shift,w)-1. );
    }
    Kn_posdt_integrand_class(RealFunction *J_,double beta_, double E_shift_, double dt_)
     : J(J_), beta(beta_), E_shift(E_shift_), dt(dt_){
    }
  };

  class Kn_negdt_integrand_class: public ComplexFunction{
  public:
    RealFunction *J;
    double beta, E_shift, dt;
    virtual std::complex<double> f(double w)const{
      if(w*w<1e-20)return 0.;
      return 2.*J->f(w)/(w*w)*(1.-cos(w*dt))*1./2.*(get_coth(beta,E_shift,w)+1. );
    }
    Kn_negdt_integrand_class(RealFunction *J_,double beta_, double E_shift_, double dt_)
     : J(J_), beta(beta_), E_shift(E_shift_), dt(dt_){
    }
  };


#endif
